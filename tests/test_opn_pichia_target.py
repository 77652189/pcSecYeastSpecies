from scripts import opn_pichia_target as opn


def test_opn_target_row_matches_project_construct():
    row = opn.load_opn_row()
    opn.validate_row(row)

    assert row["Protein name"] == "OPN"
    assert row["abbreviation"] == "OPN"
    assert row["sequence"].endswith(opn.MATURE_OPN)
    assert row["sequence"].startswith(opn.PROJECT_ALPHA_FACTOR_LEADER)
    assert int(row["Length"]) == 387


def test_opn_model_features_are_conservative_and_traceable():
    row = opn.load_opn_row()
    summary = opn.build_summary(row)

    assert summary["mature_opn_length"] == 298
    assert summary["leader_length"] == 89
    assert summary["mature_cysteine_count"] == 0
    assert summary["disulfide_sites"] == 0
    assert summary["n_glycosylation_sites_for_model"] == 0
    assert summary["o_glycosylation_sites_for_model"] == 7
    assert summary["mature_nxs_t_motifs"] == [(63, "NES"), (90, "NDS")]
    assert summary["mature_opn_internal_kex2_like_sites"] == [(159, "RR"), (231, "KR")]


def test_candidate_rows_are_model_ready_and_traceable():
    rows = opn.build_candidate_rows()
    opn.validate_candidate_rows(rows)

    by_id = {row["Protein name"]: row for row in rows}

    assert {
        "OPN_ALPHA_FULL_PROJECT",
        "OPN_ALPHA_PRE_ONLY",
        "OPN_NATIVE_SPP1",
        "OPN_OST1N23_ALPHA_PRO",
        "OPN_PPA_DDDK18",
        "OPN_PPA_PASCHR3_0030",
        "OPN_PPA_EPX1_SA",
    }.issubset(by_id)
    assert by_id["OPN_ALPHA_FULL_PROJECT"]["sequence"].startswith(opn.PROJECT_ALPHA_FACTOR_LEADER)
    assert by_id["OPN_NATIVE_SPP1"]["sequence"].startswith(opn.HUMAN_SPP1_NATIVE_SP)
    assert by_id["OPN_PPA_DDDK18"]["sp sequence"] == opn.PPA_DDDK18_SP
    assert all(str(row["sequence"]).endswith(opn.MATURE_OPN) for row in rows)


def test_candidate_metadata_flags_kex2_risk_for_all_candidates():
    metadata = opn.build_candidate_metadata()

    assert len(metadata) == len(opn.candidate_definitions())
    assert all(item["mature_opn_internal_kex2_like_sites"] for item in metadata)
    assert any(item["category"] == "pichia_native_signal" for item in metadata)


def test_candidate_csv_on_disk_matches_generator():
    disk_rows = opn.load_candidate_rows()
    generated_rows = opn.build_candidate_rows()

    assert [row["Protein name"] for row in disk_rows] == [row["Protein name"] for row in generated_rows]
    assert [row["sequence"] for row in disk_rows] == [row["sequence"] for row in generated_rows]
    assert [int(row["Length"]) for row in disk_rows] == [row["Length"] for row in generated_rows]
