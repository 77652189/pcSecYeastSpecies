# 在线外部数据库证据层架构

状态：active
最后更新：2026-07-09

## 实现状态

Round 1-9 已完成并推送：external reference schema/cache、query builder、受控在线 clients、name resolution、KO/OE gene function evidence、external GPR candidate evidence、service/Streamlit 展示、screen/report/LLM fact pack 透传和端到端 smoke 均已落地。

外部 GEM/GPR 资源优先级审计与导入主线 Round A-G 也已完成并推送：已覆盖 external model inventory、受控 artifact cache/manual-required 记录、SBML reaction/GPR association parse、GPR source priority/conflict report、当前 Pichia GEM reaction/gene mapping status，以及 service/Streamlit/report/fact pack 透传。最近本地 smoke 写入 ignored `local_runs/external_model_gpr_inventory/round_g_smoke_20260709/`。

联网 smoke、外部模型审计和映射产物仍保留在 ignored `local_runs/`，不作为稳定科学资产提交；需要人工复核后再决定是否提升 production cache。

## 目标

用受控联网方式获取外部数据库和外部模型证据，服务三类研发问题：

1. Pichia 内部 `gene_id`、locus、common name、alias 的标准化命名校对。
2. KO/OE 候选的外部 gene function、EC、GO、pathway、orthology 补充。
3. 酿酒酵母同源规则和外部 GPR candidate 的可迁移性复核。

该层只提供证据和候选判断，不直接改变 KO/OE `recommendation_tier`，不自动写回 curated catalog，不替代当前 Pichia GEM 的 GPR 可执行性或 phenotype evidence。

## 边界

- 联网只能发生在明确的 builder、refresh job 或用户手动刷新中。
- Streamlit 页面加载、筛选和导出只读 cache，不默认联网。
- 核心仿真、KO/OE screen、report generation 不依赖实时网络。
- 外部 GPR 规则只能先作为 `external_gpr_candidate`；只有映射到当前 Pichia GEM 的 gene/reaction 后，才可标记为 `model_gpr_executable`。
- 外部 GEM/GPR source priority 是证据优先级，不是规则合并器；冲突来源必须进入 manual review，不自动写回当前 Pichia GEM。
- 所有记录必须保留 `source_database`、`source_version`、`source_url`、`source_query`、`retrieved_at`、warning 和原始记录 hash。
- 初始产物写入 ignored `local_runs/`；人工复核后再讨论是否提升为稳定科学资产。

## 架构

```text
Online sources
  UniProt / NCBI / SGD
  Optional future source: KEGG
  yeast-GEM / BiGG / BioModels / ModelSEED / curated SBML
    -> external_refs clients
      -> retry / rate limit / provenance
        -> external reference cache
          -> name resolution
          -> gene function evidence
          -> GPR candidate evidence
            -> app service facade
              -> Streamlit status / manual refresh / export
```

模块职责：

- `python_pichia/src/pcsec_pichia/external_refs/`：schema、query、fetch、cache、merge、classification。
- `scripts/`：批量构建和刷新 cache，可联网。
- `app/services/`：读取 cache、提交刷新任务、返回状态和导出数据。
- `app/ui/`：展示状态、证据、冲突、人工复核理由和下载入口。

## 模块布局

```text
python_pichia/src/pcsec_pichia/external_refs/
  __init__.py
  schema.py
  cache_io.py
  queries.py
  clients.py
  uniprot.py
  sgd.py
  ncbi.py
  merge.py
  name_resolution.py
  gene_function.py
  gpr_sources.py
  gpr_candidates.py
  refresh.py

scripts/
  build_pichia_external_reference_cache.py

app/services/
  pichia_external_reference_service.py
```

## 依赖

本主线一次性声明完整依赖，避免后续轮次反复修改依赖文件：

```text
httpx
tenacity
biopython
python-libsbml
pydantic
```

职责：

- `httpx`：所有在线数据库客户端统一使用的 HTTP 层。
- `tenacity`：retry/backoff 和可测试的瞬时失败恢复。
- `biopython`：NCBI Entrez、GenBank/FASTA/sequence record 解析。
- `python-libsbml`：外部 SBML/GEM 中 reaction、gene product 和 GPR 解析。
- `pydantic`：manifest、service payload、cache record 的边界校验；内部核心对象仍可用 dataclass 表达。

暂不引入：

```text
highspy
memote
pyhmmer
bioservices
requests-cache
```

## 数据结构

```python
@dataclass(frozen=True)
class ExternalReferenceQuery:
    internal_gene_id: str
    organism: str
    taxon_id: str | None
    query_name: str
    query_aliases: tuple[str, ...]
    sequence_accession: str | None = None
    protein_sequence_sha256: str | None = None
    source_context: str = "homology_cache"
```

```python
@dataclass(frozen=True)
class ExternalFetchConfig:
    sources: tuple[str, ...] = ("uniprot", "sgd", "ncbi")
    timeout_seconds: float = 20.0
    retry_attempts: int = 3
    min_interval_seconds: float = 0.2
    user_agent: str = "pcSecPichia external reference cache builder"
    offline_cache_dir: Path | None = None
    ncbi_api_key_env: str = "NCBI_API_KEY"
```

```python
@dataclass(frozen=True)
class ExternalReferenceRecord:
    source_database: str
    source_version: str
    source_url: str
    source_query: str
    organism: str
    taxon_id: str | None
    accession: str
    gene_name: str | None
    locus_tag: str | None
    ordered_locus_name: str | None
    aliases: tuple[str, ...]
    protein_names: tuple[str, ...]
    cross_references: Mapping[str, tuple[str, ...]]
    sequence_accession: str | None
    retrieved_at: str
    http_status: int | None
    raw_record_sha256: str | None
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class ExternalGeneFunctionEvidence:
    internal_gene_id: str
    source_database: str
    accession: str
    organism: str
    gene_name: str | None
    aliases: tuple[str, ...]
    protein_names: tuple[str, ...]
    ec_numbers: tuple[str, ...]
    go_terms: tuple[str, ...]
    pathway_ids: tuple[str, ...]
    orthology_ids: tuple[str, ...]
    reviewed_status: str | None
    evidence_confidence: str
    retrieved_at: str
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class ExternalReactionAssociation:
    source_database: str
    source_model_id: str | None
    source_reaction_id: str
    source_reaction_name: str | None
    source_gene_rule: str | None
    source_gene_ids: tuple[str, ...]
    ec_numbers: tuple[str, ...]
    cross_references: Mapping[str, tuple[str, ...]]
    equation: str | None
    compartment: str | None
    retrieved_at: str
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class ExternalGprCandidateEvidence:
    pichia_gene_id: str
    query_gene_id: str
    source_database: str
    source_model_id: str | None
    source_reaction_id: str
    source_gene_rule: str | None
    mapped_model_reaction_id: str | None
    gene_mapping_status: str
    reaction_mapping_status: str
    gpr_transfer_status: str
    confidence: str
    supporting_gene_evidence: tuple[str, ...]
    blocking_reasons: tuple[str, ...]
    warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class KoOeExternalGeneEvidence:
    pichia_gene_id: str
    standard_name: str | None
    external_name_status: str
    function_evidence: tuple[ExternalGeneFunctionEvidence, ...]
    gpr_candidates: tuple[ExternalGprCandidateEvidence, ...]
    model_executable_gene_id: str | None
    model_gpr_executable: bool
    manual_review_reasons: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class ExternalReferenceCacheManifest:
    generated_at: str
    cache_schema_version: str
    query_count: int
    record_count: int
    failed_query_count: int
    source_counts: Mapping[str, int]
    input_cache_fingerprint: str | None
    warnings: tuple[str, ...]
```

## 函数契约

### schema.py / cache_io.py

```python
def write_external_reference_cache(
    records: Iterable[ExternalReferenceRecord],
    output_path: Path,
    *,
    manifest: ExternalReferenceCacheManifest,
) -> Path: ...
```

```python
def load_external_reference_cache(path: Path) -> tuple[ExternalReferenceRecord, ...]: ...
```

```python
def validate_external_reference_cache(path: Path) -> ExternalReferenceCacheManifest: ...
```

### queries.py

```python
def build_external_reference_queries_from_homology_cache(
    homology_cache_path: Path,
    *,
    include_model_genes: bool = True,
    include_sce_queries: bool = True,
) -> tuple[ExternalReferenceQuery, ...]: ...
```

```python
def build_external_reference_queries_from_gene_catalog(
    catalog_rows: Iterable[Mapping[str, Any]],
) -> tuple[ExternalReferenceQuery, ...]: ...
```

```python
def normalize_external_query_name(name: str) -> str: ...
```

### clients.py / source clients

```python
class ExternalReferenceClient(Protocol):
    source_database: str

    def fetch(self, query: ExternalReferenceQuery, config: ExternalFetchConfig) -> ExternalFetchResult: ...
```

```python
def fetch_external_references(
    queries: Iterable[ExternalReferenceQuery],
    clients: Iterable[ExternalReferenceClient],
    config: ExternalFetchConfig,
) -> tuple[ExternalFetchResult, ...]: ...
```

```python
def fetch_uniprot_reference(query: ExternalReferenceQuery, config: ExternalFetchConfig) -> ExternalFetchResult: ...
def fetch_sgd_reference(query: ExternalReferenceQuery, config: ExternalFetchConfig) -> ExternalFetchResult: ...
def fetch_ncbi_gene_reference(query: ExternalReferenceQuery, config: ExternalFetchConfig) -> ExternalFetchResult: ...
```

KEGG 不属于当前已实现 client 集；如需新增来源，应单独评估依赖、API 边界和验收命令。

### merge.py / name_resolution.py

```python
def merge_external_fetch_results(
    results: Iterable[ExternalFetchResult],
) -> tuple[ExternalReferenceRecord, ...]: ...
```

```python
def classify_external_name_consistency(
    *,
    internal_gene_id: str,
    internal_common_name: str | None,
    internal_aliases: Iterable[str],
    external_records: Iterable[ExternalReferenceRecord],
) -> NameResolutionCandidate: ...
```

```python
def attach_external_references_to_name_audit(
    name_audit_rows: Iterable[Mapping[str, Any]],
    records: Iterable[ExternalReferenceRecord],
) -> tuple[Mapping[str, Any], ...]: ...
```

### gene_function.py

```python
def build_gene_function_evidence(
    *,
    internal_gene_id: str,
    external_records: Iterable[ExternalReferenceRecord],
) -> tuple[ExternalGeneFunctionEvidence, ...]: ...
```

```python
def classify_gene_function_confidence(
    evidence: ExternalGeneFunctionEvidence,
) -> tuple[str, tuple[str, ...]]: ...
```

### gpr_sources.py / gpr_candidates.py

```python
def fetch_external_model_reaction_associations(
    *,
    source_database: str,
    model_id: str,
    gene_or_reaction_query: str,
    config: ExternalFetchConfig,
) -> tuple[ExternalReactionAssociation, ...]: ...
```

```python
def parse_sbml_gpr_associations(
    sbml_path: Path,
    *,
    source_database: str,
    source_model_id: str,
) -> tuple[ExternalReactionAssociation, ...]: ...
```

```python
def build_external_gpr_candidates(
    *,
    pichia_gene_id: str,
    query_gene_id: str,
    gene_function_evidence: Iterable[ExternalGeneFunctionEvidence],
    reaction_associations: Iterable[ExternalReactionAssociation],
    current_model_reaction_ids: Iterable[str],
    reaction_crosswalk: Mapping[str, str] | None = None,
) -> tuple[ExternalGprCandidateEvidence, ...]: ...
```

```python
def classify_gpr_transfer_status(
    *,
    gene_mapping_status: str,
    reaction_mapping_status: str,
    source_gene_rule: str | None,
    mapped_model_reaction_id: str | None,
    in_current_model_gene_index: bool,
) -> tuple[str, tuple[str, ...]]: ...
```

推荐状态：

```text
model_gpr_confirmed
external_gpr_candidate
reaction_mapping_required
gene_mapping_required
source_rule_missing
conflicting_gpr_sources
not_in_current_model
manual_review_required
```

### refresh.py

```python
def build_external_reference_cache(
    queries: Iterable[ExternalReferenceQuery],
    output_dir: Path,
    *,
    config: ExternalFetchConfig,
) -> ExternalReferenceCacheManifest: ...
```

```python
def refresh_external_reference_cache_for_homology_run(
    homology_run_dir: Path,
    *,
    output_dir: Path | None = None,
    config: ExternalFetchConfig | None = None,
) -> ExternalReferenceCacheManifest: ...
```

### app service facade

```python
def load_external_reference_status(cache_root: Path) -> dict[str, Any]: ...
def load_external_reference_browser_rows(cache_root: Path) -> list[dict[str, Any]]: ...
def submit_external_reference_refresh(*, homology_run_dir: Path, sources: tuple[str, ...], limit: int | None = None) -> dict[str, Any]: ...
def export_external_reference_rows(rows: Iterable[Mapping[str, Any]], *, file_format: str) -> bytes: ...
```

## 输出产物

```text
local_runs/pichia_external_reference_cache/<run_name>/
  external_reference_records.jsonl
  external_gene_function_evidence.jsonl
  external_reaction_associations.jsonl
  external_gpr_candidate_evidence.jsonl
  external_reference_manifest.json
  external_reference_summary.md
  failed_queries.jsonl
```

合并到同源审计时：

```text
local_runs/pichia_homology_cache/<run_name>/
  sce_to_pichia_name_audit.with_external.jsonl
  sce_to_pichia_name_audit.with_external.tsv
  rule_transfer_external_evidence.jsonl
  ko_oe_external_gene_evidence.jsonl
```

## 验收标准

默认测试不得依赖实时网络：

```powershell
python -m pytest -q python_pichia\tests\test_external_refs_schema.py python_pichia\tests\test_external_refs_queries.py python_pichia\tests\test_external_refs_clients.py
python -m pytest -q python_pichia\tests\test_external_refs_merge.py python_pichia\tests\test_external_refs_name_resolution.py
python -m pytest -q python_pichia\tests\test_external_refs_gene_function.py python_pichia\tests\test_external_refs_gpr_candidates.py
python -m pytest -q tests\test_pichia_external_reference_service_contract.py tests\test_pichia_homology_audit_service_contract.py
python -m compileall -q app python_pichia\src scripts
git diff --name-only -- Code Model Enzymedata Results requirements.txt python_pichia\pyproject.toml
```

联网 smoke 单独运行：

```powershell
python scripts\build_pichia_external_reference_cache.py --sources uniprot,sgd --limit 10 --output-dir local_runs\pichia_external_reference_cache\smoke
```

## 非目标

- 不把实时联网查询放入核心仿真路径。
- 不让外部数据库名称覆盖模型内 `gene_id`。
- 不把外部 GEM 的 GPR 规则原样写入当前 Pichia GEM。
- 不把 EC、GO、pathway 或外部 reaction association 当作当前模型已可执行 GPR。
- 不把同源证据或外部注释自动写入 curated catalog。
- 不把外部数据库注释当作 phenotype evidence。
