from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python_pichia" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from pcsec_pichia.external_refs import (
    ExternalModelArtifactRequest,
    build_artifact_requests_from_inventory,
    cache_external_model_artifacts,
    default_external_model_inventory_records,
)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    output_dir = _resolve_output_dir(args.output_dir)
    requests = list(build_artifact_requests_from_inventory(default_external_model_inventory_records()))
    requests.extend(_explicit_download_requests(args.download))
    outputs = cache_external_model_artifacts(
        tuple(requests),
        output_dir,
        timeout_seconds=args.timeout,
    )
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "manifest_path": str(outputs.manifest_path),
                "failures_path": str(outputs.failures_path),
                "request_count": outputs.manifest.request_count,
                "downloaded_count": outputs.manifest.downloaded_count,
                "failed_count": outputs.manifest.failed_count,
                "manual_required_count": outputs.manifest.manual_required_count,
                "checksum_mismatch_count": outputs.manifest.checksum_mismatch_count,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cache external GEM model artifacts into local_runs.")
    parser.add_argument(
        "--output-dir",
        default="",
        help="Output directory. Defaults to local_runs/external_model_gpr_inventory/<timestamp>/artifacts_cache.",
    )
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument(
        "--download",
        action="append",
        default=[],
        metavar="MODEL_ID=URL=FILENAME=TYPE[=SHA256]",
        help="Explicit direct artifact download. Repository or publication landing pages should not be passed here.",
    )
    return parser.parse_args(argv)


def _explicit_download_requests(values: list[str]) -> tuple[ExternalModelArtifactRequest, ...]:
    requests: list[ExternalModelArtifactRequest] = []
    for value in values:
        parts = value.split("=")
        if len(parts) not in {4, 5}:
            raise ValueError("--download must be MODEL_ID=URL=FILENAME=TYPE[=SHA256].")
        model_id, url, filename, artifact_type = parts[:4]
        expected_sha256 = parts[4] if len(parts) == 5 else ""
        requests.append(
            ExternalModelArtifactRequest(
                model_id=model_id,
                artifact_url=url,
                artifact_type=artifact_type,
                filename=filename,
                expected_sha256=expected_sha256,
            )
        )
    return tuple(requests)


def _resolve_output_dir(value: str) -> Path:
    if value:
        path = Path(value)
        return path if path.is_absolute() else REPO_ROOT / path
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return REPO_ROOT / "local_runs" / "external_model_gpr_inventory" / stamp / "artifacts_cache"


if __name__ == "__main__":
    raise SystemExit(main())
