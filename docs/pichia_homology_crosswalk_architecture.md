# BLAST/RBH 同源映射架构

状态：active  
最后更新：2026-07-07

## 目标

建立一个离线、可审计的同源映射证据层，把酿酒酵母分泌工程知识迁移到 Pichia 模型前先转化为结构化 crosswalk/cache。

该层只回答：

```text
这个 S. cerevisiae 基因在 Pichia 中是否有可信同源候选？
这个候选是否存在于当前 Pichia GEM gene_index？
是否需要人工复核？
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
  -> JSONL / TSV / Markdown report
```

## 建议模块布局

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

scripts/
  build_pichia_homology_cache.py

python_pichia/tests/
  test_homology_sequence_sources.py
  test_homology_blast_parser.py
  test_homology_rbh.py
  test_homology_crosswalk.py
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
    sce_symbol: str
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

### crosswalk.py

```python
def build_homology_crosswalk(
    sce_queries: Iterable[CatalogHomologyQuery],
    sce_records: Iterable[ProteinRecord],
    pichia_records: Iterable[ProteinRecord],
    model_gene_index: set[str],
    blast_config: BlastConfig,
) -> HomologyCrosswalk:
    """Build the full offline crosswalk from local sequence assets and BLAST results."""
```

```python
def write_homology_cache(crosswalk: HomologyCrosswalk, jsonl_path: Path, tsv_path: Path) -> CacheWriteResult:
    """Write deterministic cache outputs for review and downstream use."""
```

```python
def load_homology_cache(path: Path) -> HomologyCrosswalk:
    """Load a previously generated cache without rerunning BLAST."""
```

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
homology_cache_summary.md
blast_forward.tsv
blast_reverse.tsv
```

These outputs are review artifacts. They should not be committed until the project decides to promote a curated cache to stable scientific input.

## Integration path

Phase 1: build cache and report only.  
Phase 2: join cache into gene evidence/capability profiles as `homology_evidence`.  
Phase 3: expose homology evidence in KO/OE recommendation explanations.  
Phase 4: after manual review, promote selected mappings into curated catalog or fixture data.

At every phase, homology evidence remains separate from phenotype evidence and model executability.
