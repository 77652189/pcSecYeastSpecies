from __future__ import annotations

from pcsec_pichia.adapters.lp_parser import parse_lp_file
from pcsec_pichia.adapters.lp_writer import write_stoichiometric_lp
from pcsec_pichia.adapters.mat_loader import MatStructLoader
from pcsec_pichia.core.pichia_media import apply_glucose_reference_conditions, set_fixed_growth_rate
from pcsec_pichia.core.paths import ProjectPaths


def test_lp_writer_exports_sv_constraints_and_bounds_with_matlab_variable_numbering(tmp_path) -> None:
    model = MatStructLoader(ProjectPaths.discover()).load_pcsec_pichia_model()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_sv_bounds.lp"

    write_summary = write_stoichiometric_lp(fixed_mu, output, objective_reaction="Ex_glc_D")
    parsed = parse_lp_file(output)

    assert write_summary.objective_variable == "X545"
    assert write_summary.stoichiometric_constraint_count == 20195
    assert write_summary.metabolic_coupling_constraint_count == 0
    assert write_summary.secretory_coupling_constraint_count == 0
    assert write_summary.protein_mass_constraint_count == 0
    assert write_summary.mitochondrial_constraint_count == 0
    assert write_summary.proteasome_constraint_count == 0
    assert write_summary.ribosome_assembly_constraint_count == 0
    assert write_summary.bounds_count == 29026

    assert parsed.optimization_sense == "Maximize"
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 20195
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195}
    assert parsed.constraint_sense_counts == {"=": 20195}

    text = output.read_text(encoding="utf-8")
    assert "C1:" in text
    assert "0.100000 <= X486 <= 0.100000" in text
    assert "-1000.000000 <= X545 <= +infinity" in text


def test_lp_writer_exports_metabolic_coupling_constraints(tmp_path) -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    metabolic_enzymes = loader.load_metabolic_enzymedata()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_met_coupling.lp"

    write_summary = write_stoichiometric_lp(
        fixed_mu,
        output,
        objective_reaction="Ex_glc_D",
        metabolic_enzymes=metabolic_enzymes,
        mu=0.10,
    )
    parsed = parse_lp_file(output)

    assert write_summary.metabolic_coupling_constraint_count == 2732
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 22927
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195, "CM": 2732}
    assert parsed.constraint_sense_counts == {"=": 22927}

    text = output.read_text(encoding="utf-8")
    assert "CM1:" in text
    assert "CM2732:" in text


def test_lp_writer_exports_secretory_coupling_constraints(tmp_path) -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    metabolic_enzymes = loader.load_metabolic_enzymedata()
    secretory_enzymes = loader.load_secretory_enzymedata()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_sec_coupling.lp"

    write_summary = write_stoichiometric_lp(
        fixed_mu,
        output,
        objective_reaction="Ex_glc_D",
        metabolic_enzymes=metabolic_enzymes,
        secretory_enzymes=secretory_enzymes,
        mu=0.10,
    )
    parsed = parse_lp_file(output)

    assert write_summary.metabolic_coupling_constraint_count == 2732
    assert write_summary.secretory_coupling_constraint_count == 58
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 22985
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195, "CM": 2790}
    assert parsed.constraint_sense_counts == {"=": 22985}

    text = output.read_text(encoding="utf-8")
    assert "CM2733:" in text
    assert "CM2790:" in text


def test_lp_writer_exports_total_enzyme_and_unmodeled_er_constraints(tmp_path) -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    metabolic_enzymes = loader.load_metabolic_enzymedata()
    secretory_enzymes = loader.load_secretory_enzymedata()
    combined_enzymes = loader.load_combined_enzymedata()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_protein_mass.lp"

    write_summary = write_stoichiometric_lp(
        fixed_mu,
        output,
        objective_reaction="Ex_glc_D",
        metabolic_enzymes=metabolic_enzymes,
        secretory_enzymes=secretory_enzymes,
        combined_enzymes=combined_enzymes,
        mu=0.10,
    )
    parsed = parse_lp_file(output)

    assert write_summary.protein_mass_constraint_count == 2
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 22987
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195, "CM": 2792}
    assert parsed.constraint_sense_counts == {"=": 22987}

    by_name = {constraint.name: constraint for constraint in parsed.constraints}
    assert by_name["CM2791"].variable_count == 4212
    assert by_name["CM2792"].variable_count == 1

    text = output.read_text(encoding="utf-8")
    assert "CM2791:" in text
    assert " = 0.032552119000000" in text
    assert "CM2792: 40.000000000000000 X28977 = 0.001480000000000" in text


def test_lp_writer_exports_mitochondrial_constraint(tmp_path) -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    metabolic_enzymes = loader.load_metabolic_enzymedata()
    secretory_enzymes = loader.load_secretory_enzymedata()
    combined_enzymes = loader.load_combined_enzymedata()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_cmito.lp"

    write_summary = write_stoichiometric_lp(
        fixed_mu,
        output,
        objective_reaction="Ex_glc_D",
        metabolic_enzymes=metabolic_enzymes,
        secretory_enzymes=secretory_enzymes,
        combined_enzymes=combined_enzymes,
        mu=0.10,
        include_mitochondrial_constraint=True,
    )
    parsed = parse_lp_file(output)

    assert write_summary.mitochondrial_constraint_count == 1
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 22988
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195, "CM": 2792, "Cmito": 1}
    assert parsed.constraint_sense_counts == {"=": 22987, "<=": 1}

    by_name = {constraint.name: constraint for constraint in parsed.constraints}
    assert by_name["Cmito"].variable_count == 522

    text = output.read_text(encoding="utf-8")
    assert "Cmito:" in text
    assert " <= 0.005000000000000" in text


def test_lp_writer_exports_proteasome_constraint(tmp_path) -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    metabolic_enzymes = loader.load_metabolic_enzymedata()
    secretory_enzymes = loader.load_secretory_enzymedata()
    combined_enzymes = loader.load_combined_enzymedata()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_proteasome.lp"

    write_summary = write_stoichiometric_lp(
        fixed_mu,
        output,
        objective_reaction="Ex_glc_D",
        metabolic_enzymes=metabolic_enzymes,
        secretory_enzymes=secretory_enzymes,
        combined_enzymes=combined_enzymes,
        mu=0.10,
        include_mitochondrial_constraint=True,
        include_proteasome_constraint=True,
    )
    parsed = parse_lp_file(output)

    assert write_summary.proteasome_constraint_count == 1
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 22989
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195, "CM": 2793, "Cmito": 1}
    assert parsed.constraint_sense_counts == {"=": 22988, "<=": 1}

    by_name = {constraint.name: constraint for constraint in parsed.constraints}
    assert by_name["CM2793"].variable_count == 1418

    text = output.read_text(encoding="utf-8")
    assert "CM2793:" in text
    assert " - 600.000000000000000 X27043 = 0" in text


def test_lp_writer_exports_ribosome_assembly_constraint_with_matlab_gap_numbering(tmp_path) -> None:
    loader = MatStructLoader(ProjectPaths.discover())
    model = loader.load_pcsec_pichia_model()
    metabolic_enzymes = loader.load_metabolic_enzymedata()
    secretory_enzymes = loader.load_secretory_enzymedata()
    combined_enzymes = loader.load_combined_enzymedata()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10_ribosome_assembly.lp"

    write_summary = write_stoichiometric_lp(
        fixed_mu,
        output,
        objective_reaction="Ex_glc_D",
        metabolic_enzymes=metabolic_enzymes,
        secretory_enzymes=secretory_enzymes,
        combined_enzymes=combined_enzymes,
        mu=0.10,
        include_mitochondrial_constraint=True,
        include_proteasome_constraint=True,
        include_ribosome_assembly_constraint=True,
    )
    parsed = parse_lp_file(output)

    assert write_summary.ribosome_assembly_constraint_count == 1
    assert parsed.objective_variable == "X545"
    assert parsed.constraint_count == 22990
    assert parsed.bounds_count == 29026
    assert parsed.distinct_variable_count == 29026
    assert parsed.max_variable_index == 29026
    assert parsed.constraint_prefix_counts == {"C": 20195, "CM": 2794, "Cmito": 1}
    assert parsed.constraint_sense_counts == {"=": 22989, "<=": 1}

    by_name = {constraint.name: constraint for constraint in parsed.constraints}
    assert by_name["CM2796"].variable_count == 2

    text = output.read_text(encoding="utf-8")
    assert "CM2796: X27039 - 1200000.000000000000000 X27041 = 0" in text




def test_lp_writer_wraps_long_lines_for_soplex_compatibility(tmp_path) -> None:
    paths = ProjectPaths.discover()
    model = MatStructLoader(paths).load_pcsec_pichia_model()
    glucose = apply_glucose_reference_conditions(model, media_type=4).model
    fixed_mu = set_fixed_growth_rate(glucose, 0.10)
    output = tmp_path / "python_glucose_mu0p10.lp"

    write_stoichiometric_lp(fixed_mu, output, objective_reaction="Ex_glc_D")

    longest_line = max(len(line.rstrip("\n")) for line in output.open(encoding="utf-8"))

    assert longest_line < 8190
