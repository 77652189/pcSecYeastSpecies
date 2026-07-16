from __future__ import annotations

import json

from pcsec_pichia.core.paths import ProjectPaths

from app.services.genome_wide_screen_registry import (
    RunInfo,
    latest_runs_by_group,
    list_runs,
    older_runs_by_group,
    run_group_key,
    run_scope_family,
)
from app.ui.views.genome_wide_screen import _split_result_runs


def _run(
    run_name: str,
    *,
    scope: str = "gene",
    targets: tuple[str, ...] = ("hLF",),
    status: str = "done",
    done: int = 10,
    total: int = 10,
) -> RunInfo:
    return RunInfo(
        run_name=run_name,
        status=status,
        done=done,
        total=total,
        targets=targets,
        mode="fast",
        pid=None,
        updated_at="2026-07-06T00:00:00",
        is_stale=False,
        scope=scope,
    )


def test_run_group_key_ignores_target_order() -> None:
    a = _run("a", targets=("hLF", "OPN"))
    b = _run("b", targets=("OPN", "hLF"))
    assert run_group_key(a) == run_group_key(b)


def test_run_group_key_distinguishes_different_target_sets_same_scope() -> None:
    hlf = _run("gene_hlf", scope="gene", targets=("hLF",))
    opn = _run("gene_opn", scope="gene", targets=("OPN",))
    assert run_group_key(hlf) != run_group_key(opn)


def test_latest_runs_by_group_keeps_newest_per_group_only() -> None:
    # newest-first, matching list_runs()'s ordering
    runs = [
        _run("catalog_v3", scope="catalog", targets=("hLF", "OPN")),
        _run("catalog_v2", scope="catalog", targets=("hLF", "OPN")),
        _run("catalog_v1", scope="catalog", targets=("hLF", "OPN")),
    ]

    latest = latest_runs_by_group(runs)

    assert [run.run_name for run in latest] == ["catalog_v3"]


def test_latest_runs_by_group_does_not_collapse_different_targets() -> None:
    """Gene-scope hLF and gene-scope OPN are different analyses, not repeats of each
    other - both must survive even though they share a scope."""
    runs = [
        _run("overnight_hLF_full", scope="gene", targets=("hLF",)),
        _run("overnight_OPN_full", scope="gene", targets=("OPN",)),
    ]

    latest = latest_runs_by_group(runs)

    assert {run.run_name for run in latest} == {"overnight_hLF_full", "overnight_OPN_full"}


def test_latest_runs_by_group_does_not_let_gene_smoke_supersede_full_target_run() -> None:
    """A tiny gene-scope smoke test for hLF is not a replacement for the full 1025-gene hLF run."""
    runs = [
        _run("phase5_solver_retry_smoke", scope="gene", targets=("hLF",), done=2, total=2),
        _run("overnight_hLF_full", scope="gene", targets=("hLF",), done=1025, total=1025),
    ]

    latest = latest_runs_by_group(runs)
    older = older_runs_by_group(runs)

    assert {run.run_name for run in latest} == {"phase5_solver_retry_smoke", "overnight_hLF_full"}
    assert older == []


def test_run_scope_family_distinguishes_full_gene_sweeps_from_smoke_runs() -> None:
    full = _run("overnight_hLF_full", scope="gene", targets=("hLF",), done=1025, total=1025)
    smoke = _run("phase5_solver_retry_smoke", scope="gene", targets=("hLF",), done=2, total=2)
    catalog = _run("catalog_reaction_screen", scope="catalog", targets=("hLF", "OPN"), done=244, total=244)

    assert run_scope_family(full) == "gene"
    assert run_scope_family(smoke) == "gene_limited"
    assert run_scope_family(catalog) == "catalog"


def test_latest_runs_by_group_surfaces_in_progress_rerun_over_stale_done_one() -> None:
    """An in-progress re-run is newer (by construction, list_runs() sorts newest-first via
    mtime) than an older completed run in the same group, so it should be what shows -
    a fresh attempt superseding a stale success, not the other way around."""
    runs = [
        _run("catalog_rerun", scope="catalog", targets=("hLF", "OPN"), status="running"),
        _run("catalog_old", scope="catalog", targets=("hLF", "OPN"), status="done"),
    ]

    latest = latest_runs_by_group(runs)

    assert [run.run_name for run in latest] == ["catalog_rerun"]


def test_older_runs_by_group_is_the_complement_of_latest() -> None:
    runs = [
        _run("catalog_v3", scope="catalog", targets=("hLF", "OPN")),
        _run("catalog_v2", scope="catalog", targets=("hLF", "OPN")),
        _run("gene_hlf", scope="gene", targets=("hLF",)),
    ]

    older = older_runs_by_group(runs)

    assert [run.run_name for run in older] == ["catalog_v2"]


def test_older_runs_by_group_empty_when_every_group_has_one_run() -> None:
    runs = [
        _run("overnight_hLF_full", scope="gene", targets=("hLF",)),
        _run("overnight_OPN_full", scope="gene", targets=("OPN",)),
    ]

    assert older_runs_by_group(runs) == []


def test_result_run_split_does_not_show_superseded_done_run_as_latest() -> None:
    runs = [
        _run("catalog_rerun", scope="catalog", targets=("hLF", "OPN"), status="running"),
        _run("catalog_old", scope="catalog", targets=("hLF", "OPN"), status="done"),
        _run("gene_hlf", scope="gene", targets=("hLF",), status="done"),
    ]

    latest_done, older_done = _split_result_runs(runs)

    assert [run.run_name for run in latest_done] == ["gene_hlf"]
    assert [run.run_name for run in older_done] == ["catalog_old"]


def test_list_runs_reads_error_count_from_status_json(tmp_path) -> None:
    """The screen runner writes error_count whenever some tasks fail (e.g. every candidate
    for a target failed to solve, leaving a "done" run with zero result rows). list_runs()
    used to drop this field entirely, so the UI had no way to explain such a run instead of
    just crashing on an empty target list.
    """
    run_dir = tmp_path / "local_runs" / "ui_all_tasks_failed"
    run_dir.mkdir(parents=True)
    (run_dir / "status.json").write_text(
        json.dumps(
            {
                "status": "done",
                "done": 5,
                "total": 5,
                "targets": ["hLF"],
                "mode": "fast",
                "pid": None,
                "updated_at": "2026-07-15T18:55:23",
                "scope": "gene",
                "csv_path": str(run_dir / "gene_tradeoff_rows.csv"),
                "error_count": 5,
            }
        ),
        encoding="utf-8",
    )

    [run] = list_runs(ProjectPaths(repo_root=tmp_path))

    assert run.error_count == 5


def test_list_runs_defaults_error_count_to_zero_when_absent() -> None:
    """Legacy/backfilled status.json payloads never had this field; it must not crash
    or become None, since RunInfo.error_count is read as a plain int by the UI."""
    assert _run("legacy").error_count == 0
