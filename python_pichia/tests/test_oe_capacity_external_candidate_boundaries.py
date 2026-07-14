from __future__ import annotations

import ast
from pathlib import Path

import pcsec_pichia.oe_capacity as public_api
from pcsec_pichia.oe_capacity import external_candidates as compatibility_facade


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "python_pichia" / "src" / "pcsec_pichia"


def _module_ast(relative_path: str) -> ast.Module:
    return ast.parse((REPO_ROOT / relative_path).read_text(encoding="utf-8"))


def _imported_modules(relative_path: str) -> set[str]:
    modules: set[str] = set()
    for node in ast.walk(_module_ast(relative_path)):
        if isinstance(node, ast.ImportFrom) and node.module:
            modules.add(node.module)
        elif isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
    return modules


def test_external_candidates_is_a_thin_compatibility_facade() -> None:
    module = _module_ast(
        "python_pichia/src/pcsec_pichia/oe_capacity/external_candidates.py"
    )
    assert not any(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        for node in module.body
    )


def test_adr_module_ownership_boundaries_are_static() -> None:
    source_imports = _imported_modules(
        "python_pichia/src/pcsec_pichia/external_refs/capacity_sources.py"
    )
    assert not any(name.startswith("pcsec_pichia.oe_capacity") for name in source_imports)

    schema_imports = _imported_modules(
        "python_pichia/src/pcsec_pichia/oe_capacity/external_candidate_schema.py"
    )
    assert not {"csv", "json", "pathlib", "urllib.request"}.intersection(schema_imports)
    assert not any(
        name.startswith(("pcsec_pichia.loading", "pcsec_pichia.screens"))
        for name in schema_imports
    )

    io_imports = _imported_modules(
        "python_pichia/src/pcsec_pichia/oe_capacity/external_candidate_io.py"
    )
    assert not any(
        name.startswith(
            (
                "pcsec_pichia.external_refs.clients",
                "pcsec_pichia.external_refs.uniprot",
                "app",
            )
        )
        for name in io_imports
    )

    evaluation_imports = _imported_modules(
        "python_pichia/src/pcsec_pichia/oe_capacity/external_candidate_evaluation.py"
    )
    assert not any(
        name.startswith(
            (
                "pcsec_pichia.external_refs.clients",
                "pcsec_pichia.external_refs.uniprot",
                "app",
            )
        )
        for name in evaluation_imports
    )

    promotion_imports = _imported_modules(
        "python_pichia/src/pcsec_pichia/oe_capacity/external_candidate_promotion.py"
    )
    assert not any(name.startswith("app") for name in promotion_imports)

    audit_imports = _imported_modules(
        "python_pichia/src/pcsec_pichia/oe_capacity/external_candidate_audit.py"
    )
    assert not any(name.startswith("app") for name in audit_imports)


def test_cli_only_uses_public_core_audit_api() -> None:
    path = "python_pichia/tools/run_oe_capacity_external_candidate_audit.py"
    imports = _imported_modules(path)
    assert not any(name.startswith("app") for name in imports)
    assert "pcsec_pichia.oe_capacity.external_candidate_audit" in imports
    for node in ast.walk(_module_ast(path)):
        if isinstance(node, ast.ImportFrom):
            assert all(not alias.name.startswith("_") for alias in node.names)


def test_public_exports_reuse_split_module_objects() -> None:
    from pcsec_pichia.external_refs.capacity_sources import ExternalCapacitySource
    from pcsec_pichia.oe_capacity.external_candidate_evaluation import (
        build_capacity_candidate,
    )
    from pcsec_pichia.oe_capacity.external_candidate_io import (
        load_external_capacity_candidate_bundle,
    )
    from pcsec_pichia.oe_capacity.external_candidate_promotion import (
        promote_capacity_candidates,
    )
    from pcsec_pichia.oe_capacity.external_candidate_schema import (
        ExternalCapacityCandidate,
    )

    for name, expected in (
        ("ExternalCapacitySource", ExternalCapacitySource),
        ("ExternalCapacityCandidate", ExternalCapacityCandidate),
        ("build_capacity_candidate", build_capacity_candidate),
        ("load_external_capacity_candidate_bundle", load_external_capacity_candidate_bundle),
        ("promote_capacity_candidates", promote_capacity_candidates),
    ):
        assert getattr(public_api, name) is expected
        assert getattr(compatibility_facade, name) is expected


def test_external_source_validation_uses_a_neutral_shared_error() -> None:
    from pcsec_pichia.errors import OECapacityValidationError as SharedValidationError
    from pcsec_pichia.external_refs.capacity_sources import (
        ExternalCapacitySourceValidationError,
    )
    from pcsec_pichia.oe_capacity.schema import OECapacityValidationError

    assert ExternalCapacitySourceValidationError is SharedValidationError
    assert OECapacityValidationError is SharedValidationError


def test_audit_orchestration_is_public_and_app_independent() -> None:
    from pcsec_pichia.oe_capacity.external_candidate_audit import (
        ExternalCapacityAuditRequest,
        prepare_external_candidate_runtime,
        run_external_capacity_candidate_audit,
    )

    assert callable(prepare_external_candidate_runtime)
    assert callable(run_external_capacity_candidate_audit)
    assert ExternalCapacityAuditRequest.__module__.endswith("external_candidate_audit")
