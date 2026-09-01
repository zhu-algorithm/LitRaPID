"""Zero-dependency local web platform for the LitRaPID display loop."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

from .display_loop import BACKENDS, design_next_round, export_panel, feedback_update
from .mrna_display import SimulationParameters, simulate_mrna_display


ROOT = Path(__file__).resolve().parents[2]
WEB_ROOT = ROOT / "web"


def _csv_roundtrip(rows: list[dict[str, Any]], operation) -> list[dict[str, str]]:
    if not rows:
        return []
    with tempfile.TemporaryDirectory(prefix="litrapid-") as tmp:
        source = Path(tmp) / "input.csv"
        output = Path(tmp) / "output.csv"
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with source.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        operation(source, output)
        with output.open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))


def api_export(payload: dict[str, Any]) -> dict[str, Any]:
    backend = str(payload.get("backend", "yeast_gpi_macrocycle"))
    if backend not in BACKENDS:
        raise ValueError(f"Unsupported backend: {backend}")
    candidates = payload.get("candidates") or []
    limit = max(0, int(payload.get("limit", 96)))
    rows = _csv_roundtrip(
        candidates,
        lambda source, output: export_panel(source, backend, output, limit),
    )
    return {"backend": backend, "count": len(rows), "rows": rows}


def api_feedback(payload: dict[str, Any]) -> dict[str, Any]:
    panel = payload.get("panel") or []
    results = payload.get("results") or []
    if not panel or not results:
        raise ValueError("panel and results are required")
    limit = max(0, int(payload.get("limit", 96)))
    with tempfile.TemporaryDirectory(prefix="litrapid-feedback-") as tmp:
        panel_path, results_path, output = (Path(tmp) / name for name in ("panel.csv", "results.csv", "ranked.csv"))
        for path, rows in ((panel_path, panel), (results_path, results)):
            fields = list(dict.fromkeys(key for row in rows for key in row))
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
        feedback_update(panel_path, results_path, output, limit)
        with output.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    return {"count": len(rows), "rows": rows}


def api_design(payload: dict[str, Any]) -> dict[str, Any]:
    backend = str(payload.get("backend", "yeast_gpi_macrocycle"))
    ranked = payload.get("ranked") or []
    if backend not in BACKENDS:
        raise ValueError(f"Unsupported backend: {backend}")
    if not ranked:
        raise ValueError("ranked candidates are required")
    limit = max(0, int(payload.get("limit", 96)))
    parents = max(1, int(payload.get("parents", 8)))
    variants = max(0, int(payload.get("variants_per_parent", 12)))
    rows = _csv_roundtrip(
        ranked,
        lambda source, output: design_next_round(source, backend, output, limit, parents, variants),
    )
    return {"backend": backend, "count": len(rows), "rows": rows}


def api_simulate_mrna(payload: dict[str, Any]) -> dict[str, Any]:
    candidates = payload.get("candidates") or []
    allowed = set(SimulationParameters.__dataclass_fields__)
    supplied = {key: value for key, value in (payload.get("parameters") or {}).items() if key in allowed}
    integer_fields = {"rounds", "virtual_library_size", "ngs_reads", "seed"}
    supplied = {key: (int(value) if key in integer_fields else float(value)) for key, value in supplied.items()}
    return simulate_mrna_display(candidates, SimulationParameters(**supplied))


def api_meta() -> dict[str, Any]:
    return {
        "name": "LitRaPID in silico display",
        "version": "0.4.0",
        "backends": [
            {"id": key, "min_length": value.min_length, "max_length": value.max_length}
            for key, value in BACKENDS.items()
        ],
        "literature_route": {
            "backend": "yeast_gpi_macrocycle",
            "anchor": "cysteine-free GPI",
            "selection": "2× magnetic-bead enrichment + 4× two-color FACS",
            "readout": "Sanger + NGS",
        },
        "mrna_display_model": {
            "steps": ["puromycin fusion", "translation", "RaPID-like cyclization", "binding/wash", "RT-PCR", "NGS"],
            "default_rounds": 6,
            "virtual_library_size": 10**12,
        },
    }


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(WEB_ROOT), **kwargs)

    def _json(self, status: int, data: dict[str, Any]) -> None:
        body = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/api/meta":
            self._json(200, api_meta())
            return
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        try:
            length = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(length) or b"{}")
            routes = {"/api/export": api_export, "/api/feedback": api_feedback, "/api/design": api_design,
                      "/api/simulate-mrna": api_simulate_mrna}
            if self.path not in routes:
                self._json(404, {"error": "unknown endpoint"})
                return
            self._json(200, routes[self.path](payload))
        except (ValueError, TypeError, json.JSONDecodeError) as error:
            self._json(400, {"error": str(error)})
        except Exception as error:  # pragma: no cover - final server boundary
            self._json(500, {"error": f"internal error: {error}"})

    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write("[LitRaPID] " + fmt % args + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"LitRaPID platform: http://{args.host}:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
