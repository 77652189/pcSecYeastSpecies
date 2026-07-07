from __future__ import annotations

from pcsec_pichia.analysis.shadow_lp.backends import CobraOptlangBackend, ScipyHighsBackend, SolverBackend, SolverResult
from pcsec_pichia.analysis.shadow_lp.comparison import (
    SecretionCapacityComparisonResult,
    compare_secretion_capacity,
    normalize_result_status,
)
from pcsec_pichia.analysis.shadow_lp.constraint_builders import (
    REFERENCE_LAYER_ORDER,
    build_metabolic_coupling_block,
    build_misfolding_block,
    build_mitochondrial_block,
    build_proteasome_block,
    build_protein_mass_block,
    build_ribosome_assembly_block,
    build_ribosome_translation_block,
    build_secretory_coupling_block,
    build_shadow_constraint_blocks,
)
from pcsec_pichia.analysis.shadow_lp.constraint_spec import (
    ConstraintBlock,
    ConstraintSense,
    ConstraintSpec,
    LPProblem,
    OptimizationSense,
    ShadowConstraintConfig,
)
from pcsec_pichia.analysis.shadow_lp.lp_problem import (
    AssembledLPProblem,
    ConstraintOrderEntry,
    LPAssemblyDiagnostics,
    assemble_lp_problem,
    build_shadow_ladder_lp_problems,
    lp_problem_from_model,
)
from pcsec_pichia.analysis.shadow_lp.ladder import (
    FORMAL_SHADOW_LADDER_ORDER,
    ShadowLadderLayerResult,
    ShadowLadderResult,
    run_shadow_ladder,
)
from pcsec_pichia.analysis.shadow_lp.model_adapter import (
    ShadowTargetPreparation,
    fixed_growth_bounds,
    prepare_builtin_shadow_target,
    prepare_shadow_target,
)
from pcsec_pichia.analysis.shadow_lp.reports import (
    NO_ABSOLUTE_YIELD_STATEMENT,
    ShadowHardcodeAuditResult,
    render_shadow_ladder_markdown,
    render_shadow_ladder_report_payload,
    run_shadow_hardcode_audit,
    write_shadow_ladder_report,
)
from pcsec_pichia.analysis.shadow_lp.secretion_capacity import (
    SHADOW_CONSTRAINT_COUNT_KEYS,
    solve_shadow_secretion_capacity,
)
from pcsec_pichia.analysis.shadow_lp.validation import (
    ReferenceValidationResult,
    ShadowValidationMatrixCase,
    ShadowValidationMatrixResult,
    attach_reference_validation,
    run_shadow_validation_matrix,
    solve_pcsec_reference_for_validation,
    validate_shadow_ladder_against_reference,
)

__all__ = [
    "ConstraintBlock",
    "ConstraintSense",
    "ConstraintSpec",
    "AssembledLPProblem",
    "CobraOptlangBackend",
    "FORMAL_SHADOW_LADDER_ORDER",
    "LPProblem",
    "LPAssemblyDiagnostics",
    "NO_ABSOLUTE_YIELD_STATEMENT",
    "OptimizationSense",
    "REFERENCE_LAYER_ORDER",
    "ReferenceValidationResult",
    "SHADOW_CONSTRAINT_COUNT_KEYS",
    "SecretionCapacityComparisonResult",
    "ShadowConstraintConfig",
    "ShadowHardcodeAuditResult",
    "ShadowLadderLayerResult",
    "ShadowLadderResult",
    "ShadowTargetPreparation",
    "ShadowValidationMatrixCase",
    "ShadowValidationMatrixResult",
    "ScipyHighsBackend",
    "SolverBackend",
    "SolverResult",
    "ConstraintOrderEntry",
    "assemble_lp_problem",
    "attach_reference_validation",
    "build_metabolic_coupling_block",
    "build_misfolding_block",
    "build_mitochondrial_block",
    "build_proteasome_block",
    "build_protein_mass_block",
    "build_ribosome_assembly_block",
    "build_ribosome_translation_block",
    "build_secretory_coupling_block",
    "build_shadow_ladder_lp_problems",
    "build_shadow_constraint_blocks",
    "compare_secretion_capacity",
    "fixed_growth_bounds",
    "lp_problem_from_model",
    "normalize_result_status",
    "prepare_builtin_shadow_target",
    "prepare_shadow_target",
    "render_shadow_ladder_markdown",
    "render_shadow_ladder_report_payload",
    "run_shadow_hardcode_audit",
    "run_shadow_ladder",
    "run_shadow_validation_matrix",
    "solve_pcsec_reference_for_validation",
    "solve_shadow_secretion_capacity",
    "validate_shadow_ladder_against_reference",
    "write_shadow_ladder_report",
]
