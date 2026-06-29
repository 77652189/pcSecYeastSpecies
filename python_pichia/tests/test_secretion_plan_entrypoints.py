from __future__ import annotations

from pathlib import Path

from pcsec_pichia.secretion_plan import (
    SecretionPlanResult,
    build_secretion_plan,
    is_supported_target,
    summarize_secretion_plan,
)
from pcsec_pichia.targets import TargetSpec, load_builtin_targets, load_custom_targets_json


REPO_ROOT = Path(__file__).resolve().parents[2]


def _by_id(targets: list[TargetSpec]) -> dict[str, TargetSpec]:
    return {target.target_id: target for target in targets}


def test_opn_builtin_builds_og_secretion_plan() -> None:
    opn = _by_id(load_builtin_targets(REPO_ROOT))["OPN_ALPHA_FULL_PROJECT"]

    result = build_secretion_plan(opn)

    assert isinstance(result, SecretionPlanResult)
    assert result.supported is True
    assert result.route_kind == "opn_like_soluble_secretory"
    assert result.ptm_counts["o_glycosylation_sites"] == 7
    assert result.ptm_counts["disulfide_sites"] == 0
    assert result.ptm_counts["n_glycosylation_sites"] == 0
    assert result.stage_counts["folding"] >= 1
    assert result.stage_counts["golgi_processing"] >= 1
    assert any("_OG_EROG_" in reaction_id for reaction_id in result.reaction_ids)
    assert any("_GLOG_Golgi_O_linked_manosylation" in reaction_id for reaction_id in result.reaction_ids)
    assert result.raw_plan["reaction_count"] == result.reaction_count


def test_opn_leader_candidates_build_og_secretion_plans() -> None:
    targets = _by_id(load_builtin_targets(REPO_ROOT))

    for target_id in ("OPN_ALPHA_FULL_PROJECT", "OPN_PPA_DDDK18", "OPN_PPA_PASCHR3_0030"):
        result = build_secretion_plan(targets[target_id])

        assert result.supported is True
        assert result.route_kind == "opn_like_soluble_secretory"
        assert result.ptm_counts["o_glycosylation_sites"] == 7
        assert any("_OG_EROG_" in reaction_id for reaction_id in result.reaction_ids)
        assert any("_GLOG_Golgi_O_linked_manosylation" in reaction_id for reaction_id in result.reaction_ids)
        assert f"{target_id} exchange" in result.reaction_ids


def test_hlf_builtin_builds_soluble_secretory_plan_with_dsb_and_ng() -> None:
    hlf = _by_id(load_builtin_targets(REPO_ROOT))["hLF"]

    result = build_secretion_plan(hlf)
    summary = summarize_secretion_plan(hlf)

    assert result.supported is True
    assert result.route_kind == "soluble_secretory"
    assert result.ptm_counts["disulfide_sites"] == 21
    assert result.ptm_counts["n_glycosylation_sites"] == 4
    assert result.ptm_counts["o_glycosylation_sites"] == 0
    assert any("_DSB_" in reaction_id for reaction_id in result.reaction_ids)
    assert any("_ERNG_NG_" in reaction_id for reaction_id in result.reaction_ids)
    assert any("_GLNG_Golgi_N_linked_glycosylation" in reaction_id for reaction_id in result.reaction_ids)
    assert "hLF exchange" in result.reaction_ids
    assert summary["target_id"] == "hLF"
    assert summary["route_kind"] == "soluble_secretory"
    assert summary["ptm_counts"] == result.ptm_counts


def test_custom_json_targets_build_matching_secretion_plans() -> None:
    targets = _by_id(load_custom_targets_json(REPO_ROOT / "local_runs" / "pichia_hlf_opn_probe" / "targets.example.json"))

    opn_custom = build_secretion_plan(targets["OPN_CUSTOM"])
    hlf_custom = build_secretion_plan(targets["HLF_CUSTOM"])

    assert opn_custom.supported is True
    assert opn_custom.route_kind == "opn_like_soluble_secretory"
    assert opn_custom.ptm_counts["o_glycosylation_sites"] == 7
    assert any("_OG_EROG_" in reaction_id for reaction_id in opn_custom.reaction_ids)

    assert hlf_custom.supported is True
    assert hlf_custom.route_kind == "soluble_secretory"
    assert hlf_custom.ptm_counts["disulfide_sites"] == 21
    assert hlf_custom.ptm_counts["n_glycosylation_sites"] == 4
    assert any("_DSB_" in reaction_id for reaction_id in hlf_custom.reaction_ids)
    assert any("_ERNG_NG_" in reaction_id for reaction_id in hlf_custom.reaction_ids)


def test_unsupported_target_is_explicitly_marked() -> None:
    target = TargetSpec(
        target_id="TM_TEST",
        protein_id="TM_TEST",
        mature_sequence="ACDEFGHIK",
        leader_sequence="MMAA",
        signal_peptide_sequence="MM",
        through_er=True,
        localization="e",
        disulfide_sites=0,
        n_glycosylation_sites=0,
        o_glycosylation_sites=0,
        transmembrane=1,
        gpi_sites=0,
        cotranslation=0,
        source="test",
    )

    result = build_secretion_plan(target)

    assert is_supported_target(target) is False
    assert result.supported is False
    assert result.route_kind == "unsupported"
    assert result.ptm_counts["transmembrane"] == 1
