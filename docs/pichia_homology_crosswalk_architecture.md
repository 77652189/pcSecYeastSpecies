# BLAST/RBH 同源映射架构

状态：active  
最后更新：2026-07-07

## 目标

建立一个离线、可审计的同源映射证据层，把酿酒酵母分泌工程知识迁移到 Pichia 模型前先转化为结构化 crosswalk/cache，并把这些结果作为 Streamlit 中“基因命名标准化 + 同源规则迁移评估”的只读产品入口。

该层只回答：

```text
这个 S. cerevisiae 基因在 Pichia 中是否有可信同源候选？
这个候选是否存在于当前 Pichia GEM gene_index？
是否需要人工复核？
内部 common_name / gene_id 和外部/跨物种名称是否一致？
同源规则迁移是否可作为解释层进入人工复核？
```

该层不回答：

```text
敲除或过表达该基因是否一定提高分泌？
这个候选是否已经是 experiment_calibrated？
真实发酵产量是多少？
```

## 设计边界

- 不联网运行，不把数据库查询放入 app runtime。
- 不修改 `Code/`、`Model/`、`Enzymedata/`、`Results/`。
- 不自动修改 `SECRETION_GENE_CATALOG`。
- 不把 BLAST/RBH 命中直接当作模型可操作 gene。
- 不把 RBH、identity 或 annotation 单独升级为 phenotype evidence。
- Streamlit runtime 默认只读 cache，不默认联网、不默认运行 BLAST。
- `app/services` 只做 cache 读取、过滤、summary 和导出 facade，不承载核心科学判断。
- `app/ui` 只展示 python_pichia / service 产出的结构化结果，不实现同源或表型判断。
- 首轮 cache 输出到 `local_runs/`，验证稳定后再讨论是否升级为稳定资产。

## 数据流

```text
SECRETION_GENE_CATALOG / query list
  -> SCE symbol / alias normalization
  -> SCE ORF resolution
  -> SCE protein FASTA
  -> Pichia protein FASTA
  -> blastp SCE -> Pichia
  -> blastp Pichia -> SCE
  -> reciprocal best hit
  -> model gene_index join
  -> review_status
  -> homology cache / name audit / rule-transfer audit
  -> app service read-only facade
  -> Streamlit browser / filters / export
```

## 实际模块布局

```text
python_pichia/src/pcsec_pichia/homology/
  __init__.py
  sequence_sources.py
  blast_runner.py
  rbh.py
  crosswalk.py
  cache_schema.py
  review_rules.py
  catalog_inputs.py

python_pichia/src/pcsec_pichia/services/
  homology_evidence.py

app/services/
  pichia_homology_audit_service.py

app/ui/views/
  homology_audit.py

app/ui/
  common.py
  streamlit_app.py

scripts/
  build_pichia_homology_cache.py

python_pichia/tests/
  test_homology_sequence_sources.py
  test_homology_blast_parser.py
  test_homology_rbh.py
  test_homology_crosswalk.py
  test_homology_review_rules.py
  test_screens_entrypoints.py
  test_yield_improvement_recommendations.py

tests/
  test_pichia_homology_audit_service_contract.py
  test_pichia_homology_audit_ui_contract.py
  test_pichia_secretion_service_contract.py
```

## 数据结构

```python
@dataclass(frozen=True)
class ProteinRecord:
    organism: str
    gene_id: str
    symbol: str | None
    aliases: tuple[str, ...]
    sequence: str
    source: str
```

```python
@dataclass(frozen=True)
class BlastHit:
    query_id: str
    subject_id: str
    identity_pct: float
    alignment_length: int
    evalue: float
    bitscore: float
    query_coverage: float
    subject_coverage: float
```

```python
@dataclass(frozen=True)
class ReciprocalBestHit:
    query_id: str
    subject_id: str
    is_rbh: bool
    forward_hit: BlastHit | None
    reverse_hit: BlastHit | None
    failure_reason: str | None
```

```python
@dataclass(frozen=True)
class HomologyCrosswalkRow:
    internal_common_name: str
    query_symbol: str
    sce_orf: str | None
    pichia_gene_id: str | None
    pichia_model_gene_id: str | None
    is_rbh: bool
    identity_pct: float | None
    evalue: float | None
    query_coverage: float | None
    subject_coverage: float | None
    in_model_gene_index: bool
    review_status: str
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class NameAuditRow:
    internal_gene_id: str
    internal_common_name: str
    internal_sequence_id: str
    external_accession: str
    external_gene_name: str
    external_locus_tag: str
    external_aliases: tuple[str, ...]
    name_consistency_status: str
    review_status: str
    external_crosscheck_status: str = "not_available"
    external_crosscheck_sources: tuple[str, ...] = ()
    external_crosscheck_warnings: tuple[str, ...] = ()
```

```python
@dataclass(frozen=True)
class RuleTransferAuditRow:
    internal_common_name: str
    query_symbol: str
    sce_orf: str
    pichia_gene_id: str
    pichia_model_gene_id: str
    homology_review_status: str
    rule_transfer_status: str
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class ExternalNameReference:
    source_database: str
    source_version: str
    taxon: str
    accession: str
    gene_name: str
    locus_tag: str
    aliases: tuple[str, ...]
    retrieved_at: str
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class ExternalDatabaseCrosscheck:
    source_database: str
    source_version: str
    taxon: str
    accession: str
    gene_name: str
    locus_tag: str
    aliases: tuple[str, ...]
    retrieved_at: str
    match_status: str
    warnings: tuple[str, ...]
```

## 函数设计

### sequence_sources.py

```python
def load_protein_sequences_from_mat(path: Path, organism: str) -> tuple[ProteinRecord, ...]:
    """Load protein records from local pcSec protein sequence MAT assets."""
```

```python
def load_pichia_model_gene_index(root: Path | None = None) -> set[str]:
    """Return model gene ids currently executable through the Pichia GEM."""
```

```python
def write_fasta(records: Iterable[ProteinRecord], path: Path) -> Path:
    """Write deterministic FASTA used by local BLAST+ runs."""
```

### catalog_inputs.py

```python
def secretion_catalog_sce_queries() -> tuple[CatalogHomologyQuery, ...]:
    """Build SCE query symbols and aliases from the curated secretion catalog."""
```

```python
def normalize_sce_symbol(symbol: str) -> str:
    """Normalize common-name casing and separators without guessing biology."""
```

### blast_runner.py

```python
def find_blastp_executable(config: BlastConfig | None = None) -> Path | None:
    """Locate local blastp; return None when unavailable."""
```

```python
def make_blast_db(fasta_path: Path, db_prefix: Path, blast_bin: Path | None = None) -> BlastDbResult:
    """Create a local BLAST database from FASTA."""
```

```python
def run_blastp(query_fasta: Path, db_prefix: Path, out_tsv: Path, config: BlastConfig) -> BlastRunResult:
    """Run local blastp with deterministic tabular output."""
```

```python
def parse_blast_tsv(path: Path) -> tuple[BlastHit, ...]:
    """Parse BLAST outfmt rows into typed hit records."""
```

### rbh.py

```python
def best_hits_by_query(hits: Iterable[BlastHit]) -> dict[str, BlastHit]:
    """Pick best hit by evalue, bitscore, identity, and coverage."""
```

```python
def compute_reciprocal_best_hits(
    forward_hits: Iterable[BlastHit],
    reverse_hits: Iterable[BlastHit],
) -> tuple[ReciprocalBestHit, ...]:
    """Return RBH calls for SCE -> Pichia candidates."""
```

### review_rules.py

```python
def classify_homology_review_status(
    *,
    is_rbh: bool,
    identity_pct: float | None,
    query_coverage: float | None,
    subject_coverage: float | None,
    in_model_gene_index: bool,
    paralog_count: int = 0,
) -> tuple[str, tuple[str, ...]]:
    """Classify homolog evidence without turning it into phenotype evidence."""
```

Recommended statuses:

```text
model_ready_rbh_high_confidence
rbh_not_in_model
low_identity_review_required
coverage_review_required
paralog_risk_review_required
no_reciprocal_hit
unresolved_query_symbol
```

Additional audit statuses:

```text
name_confirmed_by_rbh
alias_confirmed_by_rbh
sequence_name_conflict
external_name_missing
internal_name_missing

rule_transfer_ready
rule_transfer_supported_not_model_operable
rule_transfer_low_confidence
rule_transfer_paralog_risk
rule_transfer_unresolved
rule_transfer_not_supported

not_available
external_match_confirmed
external_alias_confirmed
external_locus_confirmed
external_conflict
external_reference_incomplete
```

### crosswalk.py

```python
def build_homology_crosswalk(
    sce_queries: tuple[CatalogHomologyQuery, ...],
    sce_records: tuple[ProteinRecord, ...],
    pichia_records: tuple[ProteinRecord, ...],
    model_gene_index: set[str],
    rbh_calls: tuple[ReciprocalBestHit, ...],
    blast_config: BlastConfig,
) -> tuple[HomologyCrosswalkRow, ...]:
    """Build the full offline crosswalk from local sequence assets and BLAST results."""
```

```python
def build_name_audit_rows(
    crosswalk: tuple[HomologyCrosswalkRow, ...],
    external_references: tuple[ExternalNameReference, ...] = (),
) -> tuple[NameAuditRow, ...]:
    """Build name-audit rows while keeping model operability separate."""
```

```python
def build_rule_transfer_audit_rows(crosswalk: tuple[HomologyCrosswalkRow, ...]) -> tuple[RuleTransferAuditRow, ...]:
    """Build rule-transfer audit rows without changing phenotype evidence."""
```

```python
def load_external_name_reference_cache(path: Path) -> tuple[ExternalNameReference, ...]:
    """Load offline external name references from JSONL or TSV without network access."""
```

```python
def write_homology_cache(crosswalk: tuple[HomologyCrosswalkRow, ...], jsonl_path: Path, tsv_path: Path) -> CacheWriteResult:
    """Write deterministic cache outputs for review and downstream use."""
```

```python
def load_homology_cache(path: Path) -> tuple[HomologyCrosswalkRow, ...]:
    """Load a previously generated cache without rerunning BLAST."""
```

### app service facade

```python
def load_homology_audit_browser_data(...) -> dict[str, Any]:
    """Read cached homology audit rows for Streamlit display without running BLAST."""
```

```python
def export_homology_audit_rows(rows: list[dict[str, Any]], *, file_format: str = "tsv") -> bytes:
    """Export currently filtered rows as UTF-8 TSV/CSV bytes."""
```

### Streamlit UI

- `app/ui/common.py` 注册导航项：`基因命名与同源规则审计`。
- `app/ui/streamlit_app.py` 调用 `render_homology_audit()`。
- `app/ui/views/homology_audit.py` 展示 summary metrics、命名标准化、同源规则迁移评估、缓存状态与导出三个 tab。
- UI 只读 service 返回的结构化结果，不导入 `pcsec_pichia` engine，不运行 BLAST。

## Review status interpretation

| Status | Meaning | Downstream use |
| --- | --- | --- |
| `model_ready_rbh_high_confidence` | RBH passes thresholds and Pichia gene exists in model gene index | May support manual review and model execution checks |
| `rbh_not_in_model` | RBH exists but target is not model-operable | Homology evidence only |
| `low_identity_review_required` | RBH exists but identity below threshold | Manual review |
| `coverage_review_required` | Coverage too low for confident transfer | Manual review |
| `paralog_risk_review_required` | Multiple close hits or non-RBH best hit pattern | Manual review |
| `no_reciprocal_hit` | No RBH | Do not use as direct transfer evidence |
| `unresolved_query_symbol` | Catalog/common name could not resolve to SCE ORF | Fix query/crosswalk first |

## Output artifacts

Initial outputs should go under `local_runs/pichia_homology_cache/<run_name>/`:

```text
sce_to_pichia_homology_cache.jsonl
sce_to_pichia_homology_cache.tsv
sce_to_pichia_name_audit.jsonl
sce_to_pichia_name_audit.tsv
sce_to_pichia_rule_transfer_audit.jsonl
sce_to_pichia_rule_transfer_audit.tsv
homology_audit_summary.json
homology_audit_summary.md
blast_forward.tsv
blast_reverse.tsv
```

These outputs are review artifacts. They should not be committed until the project decides to promote a curated cache to stable scientific input.

Optional input for future external database review:

```text
--external-name-reference-cache <offline-jsonl-or-tsv>
```

The expected reference fields are `source_database`, `source_version`, `taxon`, `accession`, `gene_name`, `locus_tag`, `aliases`, `retrieved_at`, and `warnings`. This is an offline cache contract only; no online UniProt / NCBI / SGD / KEGG fetcher is implemented in the Streamlit runtime.

## Integration path

- Round 0: target contract and current-state audit completed.
- Round 1: homology / name / rule-transfer audit cache outputs completed.
- Round 2: app service facade completed.
- Round 3: Streamlit audit browser completed.
- Round 4: homology evidence joined into KO/OE explanations without changing recommendation-tier science boundaries.
- Round 5: offline external database crosscheck contract completed.
- Round 6: final docs and usage notes completed in active docs.

At every phase, homology evidence remains separate from phenotype evidence and model executability.

## Current usage

```powershell
python scripts\build_pichia_homology_cache.py --catalog-only --output-dir local_runs\pichia_homology_cache\manual_review
```

Then open Streamlit and use the “基因命名与同源规则审计” navigation item. The page reads the latest valid cache run, displays missing-cache guidance when cache files are absent, and exports only the currently filtered rows.
