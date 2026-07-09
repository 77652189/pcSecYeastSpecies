from __future__ import annotations

from pcsec_pichia.external_refs import (
    load_external_reference_cache,
    parse_external_gpr_artifacts,
    parse_sbml_gpr_associations,
)


def test_parse_sbml_gpr_associations_keeps_associations_candidate_only(tmp_path) -> None:
    sbml_path = _toy_sbml(tmp_path)

    associations = parse_sbml_gpr_associations(
        sbml_path,
        source_database="toy-gem",
        source_model_id="toy_model",
    )

    assert len(associations) == 1
    assert associations[0].external_reaction_id == "r_1234"
    assert associations[0].external_gene_ids == ("YBR160W", "YDR123C")
    assert associations[0].gene_rule == "(YBR160W and YDR123C)"
    assert associations[0].association_status == "external_gpr_candidate"
    assert associations[0].mapped_pichia_reaction_id is None
    assert associations[0].mapped_pichia_gene_ids == ()


def test_parse_external_gpr_artifacts_writes_association_cache_and_report(tmp_path) -> None:
    sbml_path = _toy_sbml(tmp_path)
    unsupported = tmp_path / "model.xlsx"
    unsupported.write_text("not parsed", encoding="utf-8")

    outputs = parse_external_gpr_artifacts(
        (sbml_path, unsupported),
        tmp_path / "out",
        source_database="toy-gem",
        source_model_id="toy_model",
    )

    assert outputs.association_count == 1
    assert outputs.unsupported_count == 1
    assert outputs.associations_path.name == "external_reaction_associations.jsonl"
    loaded = load_external_reference_cache(outputs.associations_path)
    assert loaded[0].record_type == "reaction_association"
    report = outputs.report_path.read_text(encoding="utf-8")
    assert "External GPR Parse Report" in report
    assert "unsupported_format" in report
    assert "model.xlsx" in report


def _toy_sbml(tmp_path):
    sbml_path = tmp_path / "toy.xml"
    sbml_path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2">
  <model id="toy_model">
    <fbc:listOfGeneProducts>
      <fbc:geneProduct fbc:id="G_YBR160W" fbc:label="YBR160W"/>
      <fbc:geneProduct fbc:id="G_YDR123C" fbc:label="YDR123C"/>
    </fbc:listOfGeneProducts>
    <listOfReactions>
      <reaction id="r_1234" name="example reaction">
        <fbc:geneProductAssociation>
          <fbc:and>
            <fbc:geneProductRef fbc:geneProduct="G_YBR160W"/>
            <fbc:geneProductRef fbc:geneProduct="G_YDR123C"/>
          </fbc:and>
        </fbc:geneProductAssociation>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
""",
        encoding="utf-8",
    )
    return sbml_path
