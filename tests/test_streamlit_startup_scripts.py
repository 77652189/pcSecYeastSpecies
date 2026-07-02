from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _read_repo_text(relative_path: str) -> str:
    return (REPO_ROOT / relative_path).read_text(encoding="utf-8")


def test_run_streamlit_defaults_to_python_draft_app_on_8502() -> None:
    script = _read_repo_text("run_streamlit.ps1")

    assert "[int]$Port = 8502" in script
    assert '"app/ui/streamlit_app.py"' in script
    assert "python_pichia\\src" in script
    assert "PYTHONPATH" in script
    assert "streamlit run $appPath" in script
    assert "8501" not in script


def test_lan_launcher_uses_current_streamlit_entrypoint_and_health_check() -> None:
    script = _read_repo_text("scripts/start_pcSecYeastSpecies_lan.ps1")

    assert "$port = 8502" in script
    assert "$appPath = \"app/ui/streamlit_app.py\"" in script
    assert "python_pichia\\src" in script
    assert "PYTHONPATH" in script
    assert "http://127.0.0.1:$port/_stcore/health" in script
    assert "app/ui/streamlit_app.py" in script
    assert "app\\ui\\streamlit_app.py" in script
    assert "& $python -m streamlit run $appPath" in script
    assert "8501" not in script


def test_desktop_shortcut_repair_points_to_8502_lan_launcher() -> None:
    script = _read_repo_text("scripts/repair_pcSecYeastSpecies_desktop_shortcut.ps1")
    launcher = _read_repo_text("start_pcSecYeastSpecies_lan.bat")

    assert '"pcSecYeastSpecies LAN 8502.lnk"' in script
    assert "port 8502" in script
    assert "start_pcSecYeastSpecies_lan.bat" in script
    assert "scripts\\start_pcSecYeastSpecies_lan.ps1" in launcher
    assert "8501" not in script


def test_streamlit_page_icon_is_repo_asset() -> None:
    source = _read_repo_text("app/ui/streamlit_app.py")
    icon_path = REPO_ROOT / "app" / "ui" / "assets" / "pcsecyeast_8502.png"

    assert icon_path.exists()
    assert icon_path.stat().st_size > 0
    assert "pcsecyeast_8502.png" in source
    assert "assets" in source
    assert "LauncherIcons" not in source
    assert "C:\\Users" not in source


def test_orphan_python_cleanup_script_is_repo_scoped_and_dry_run_by_default() -> None:
    script = _read_repo_text("scripts/stop_pcsec_orphan_python.ps1")

    assert "[switch]$Stop" in script
    assert "Dry run only" in script
    assert "python(?:\\.exe)?\"?\\s+-\\s*$" in script
    assert "WorkingDirectory" in script
    assert "Get-NetTCPConnection" in script
    assert "Stop-Process -Id $suspect.ProcessId" in script
