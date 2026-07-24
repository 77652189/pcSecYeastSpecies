from __future__ import annotations

import pytest

from pcsec_pichia.wetlab_private import (
    DEFAULT_PRIVATE_DIR_ENV,
    PrivateDataGuardError,
    private_data_available,
    read_private_text,
    resolve_private_data_dir,
    resolve_private_file,
)


def _repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    return repo


def test_unset_and_no_sibling_degrades_to_none(tmp_path) -> None:
    repo = _repo(tmp_path)
    assert resolve_private_data_dir(env={}, repo_root=repo) is None
    assert private_data_available(env={}, repo_root=repo) is False


def test_env_var_outside_repo_resolves(tmp_path) -> None:
    repo = _repo(tmp_path)
    private = tmp_path / "private"  # sibling of repo -> outside repo tree
    private.mkdir()
    env = {DEFAULT_PRIVATE_DIR_ENV: str(private)}
    assert resolve_private_data_dir(env=env, repo_root=repo) == private.resolve()


def test_env_var_inside_repo_raises_guard(tmp_path) -> None:
    repo = _repo(tmp_path)
    inside = repo / "sneaky_private"  # DANGER: private path inside the repo tree
    env = {DEFAULT_PRIVATE_DIR_ENV: str(inside)}
    with pytest.raises(PrivateDataGuardError):
        resolve_private_data_dir(env=env, repo_root=repo)


def test_conventional_sibling_fallback_used_when_unset(tmp_path) -> None:
    repo = _repo(tmp_path)
    sibling = tmp_path / "pcSec_wetlab_private"  # repo.parent / CONVENTIONAL name
    sibling.mkdir()
    assert resolve_private_data_dir(env={}, repo_root=repo) == sibling.resolve()


def test_read_private_text_graceful_when_unconfigured(tmp_path) -> None:
    repo = _repo(tmp_path)
    assert read_private_text("anything.txt", env={}, repo_root=repo) is None


def test_read_private_text_reads_file_and_missing_returns_none(tmp_path) -> None:
    repo = _repo(tmp_path)
    private = tmp_path / "private"
    private.mkdir()
    (private / "note.txt").write_text("mu-anchor-abstraction", encoding="utf-8")
    env = {DEFAULT_PRIVATE_DIR_ENV: str(private)}

    assert read_private_text("note.txt", env=env, repo_root=repo) == "mu-anchor-abstraction"
    assert read_private_text("missing.txt", env=env, repo_root=repo) is None


def test_path_traversal_out_of_private_area_is_refused(tmp_path) -> None:
    repo = _repo(tmp_path)
    private = tmp_path / "private"
    private.mkdir()
    env = {DEFAULT_PRIVATE_DIR_ENV: str(private)}
    with pytest.raises(PrivateDataGuardError):
        resolve_private_file("../escape.txt", env=env, repo_root=repo)
