from __future__ import annotations

import ast
import inspect
from pathlib import Path

from fastapi.testclient import TestClient

import app.api.pichia_secretion_api as api


REPO_ROOT = Path(__file__).resolve().parents[1]


def _imported_modules_from_source(source: str) -> set[str]:
    imported_modules: set[str] = set()
    module_ast = ast.parse(source)
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.add(node.module)
        elif isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
    return imported_modules


def test_fastapi_facade_is_experimental_and_uses_service_boundary() -> None:
    client = TestClient(api.create_app())
    openapi = client.get("/openapi.json").json()
    assert "Experimental" in openapi["info"]["title"]
    assert "service facade" in openapi["info"]["description"]

    imported_modules = _imported_modules_from_source(inspect.getsource(api))

    assert "app.services.pichia_secretion_service" in imported_modules
    assert not any(name.startswith("pcsec_pichia") for name in imported_modules)
    assert not any(name.startswith("python_pichia") for name in imported_modules)


def test_fastapi_imports_only_the_experimental_service_facade() -> None:
    module_ast = ast.parse(inspect.getsource(api))
    app_service_imports: dict[str, list[str]] = {}
    for node in ast.walk(module_ast):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("app.services"):
            app_service_imports.setdefault(node.module, []).extend(alias.name for alias in node.names)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("app.services"):
                    app_service_imports.setdefault(alias.name, []).append(alias.asname or alias.name)

    assert app_service_imports == {
        "app.services.pichia_secretion_service": [
            "SecretionRunRequest",
            "poll_background_simulation",
            "status_path_for_background_task",
            "submit_background_simulation",
        ]
    }


def test_fastapi_package_stays_at_service_facade_boundary() -> None:
    boundary_violations: list[str] = []
    allowed_app_service_modules = {"app.services.pichia_secretion_service"}
    forbidden_prefixes = (
        "pcsec_pichia",
        "python_pichia",
        "app.ui",
        "app.engines",
        "app.adapters",
    )

    for path in (REPO_ROOT / "app" / "api").rglob("*.py"):
        imported_modules = _imported_modules_from_source(path.read_text(encoding="utf-8"))
        for module_name in sorted(imported_modules):
            if module_name.startswith(forbidden_prefixes):
                boundary_violations.append(f"{path.relative_to(REPO_ROOT)}: {module_name}")
            if (
                module_name.startswith("app.services")
                and module_name not in allowed_app_service_modules
            ):
                boundary_violations.append(f"{path.relative_to(REPO_ROOT)}: {module_name}")

    assert boundary_violations == []


def test_fastapi_run_status_and_result_endpoints(monkeypatch, tmp_path: Path) -> None:
    captured = {}
    status_path = tmp_path / "status.json"

    def fake_submit(request):
        captured["request"] = request
        return "task-123", status_path

    def fake_poll(path):
        captured["poll_path"] = path
        return (
            "done",
            "仿真完成。",
            {
                "success": True,
                "target_id": "OPN_ALPHA_FULL_PROJECT",
                "result_status": "corrected_condition",
                "matlab_alignment_status": "pending",
            },
        )

    monkeypatch.setattr(api, "submit_background_simulation", fake_submit)
    monkeypatch.setattr(api, "poll_background_simulation", fake_poll)
    monkeypatch.setattr(api, "_status_path_for_task", lambda task_id: status_path)

    client = TestClient(api.create_app())

    create_response = client.post(
        "/pichia/secretion/runs",
        json={
            "target_source": "builtin",
            "target_id": "OPN_ALPHA_FULL_PROJECT",
            "enable_ribosome_translation_constraint": True,
            "enable_misfolding_constraint": True,
            "ko_gene_ids": ["PAS_chr1-4_0586"],
            "oe_reaction_ids": ["r_4041"],
            "growth_points": [0.1],
            "carbon_source_id": "methanol",
        },
    )
    assert create_response.status_code == 200
    create_payload = create_response.json()
    assert create_payload["task_id"] == "task-123"
    assert create_payload["status_url"] == "/pichia/secretion/runs/task-123/status"
    assert create_payload["result_url"] == "/pichia/secretion/runs/task-123/result"
    assert captured["request"].target_id == "OPN_ALPHA_FULL_PROJECT"
    assert captured["request"].ko_gene_ids == ("PAS_chr1-4_0586",)
    assert captured["request"].oe_reaction_ids == ("r_4041",)
    assert captured["request"].enable_misfolding_constraint is True
    assert captured["request"].carbon_source_id == "methanol"

    status_response = client.get("/pichia/secretion/runs/task-123/status")
    assert status_response.status_code == 200
    assert status_response.json() == {
        "task_id": "task-123",
        "status": "done",
        "message": "仿真完成。",
        "has_result": True,
    }

    result_response = client.get("/pichia/secretion/runs/task-123/result")
    assert result_response.status_code == 200
    result_payload = result_response.json()
    assert result_payload["status"] == "done"
    assert result_payload["result"]["target_id"] == "OPN_ALPHA_FULL_PROJECT"
    assert captured["poll_path"] == status_path
