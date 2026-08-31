"""Request / response models for the final_pipeline ranking API."""

from pydantic import BaseModel, Field


class GeneListItem(BaseModel):
    symbol: str
    ensg: str


class GeneListResponse(BaseModel):
    genes: list[GeneListItem]


class RankRequest(BaseModel):
    genes: list[str] = Field(..., min_length=1, max_length=10)
    lineage: list[str] = Field(default_factory=list, max_length=5)
    floor: float = Field(default=0.50, ge=0.0, le=1.0)
    top_n: int = Field(default=30, ge=1, le=200)


class RankLine(BaseModel):
    model_id: str
    name: str | None = None
    joint_score: float
    limiting_gene: str | None = None
    scores: dict[str, float]
    metadata: dict = {}


class RankResponse(BaseModel):
    genes: list[str]
    lineage: list[str]
    floor: float
    total_passing: int
    lines: list[RankLine]
    lineage_distribution: dict[str, int] = {}


class ExcludeLine(BaseModel):
    model_id: str
    name: str | None = None
    score_a: float
    score_b: float
    selectivity: float
    metadata: dict = {}


class ExcludeResponse(BaseModel):
    gene_a: str
    gene_b: str
    lineage: list[str]
    total_ranked: int
    lines: list[ExcludeLine]
    lineage_distribution: dict[str, int] = {}


class ExcludeManyLine(BaseModel):
    model_id: str
    name: str | None = None
    score_a: float
    exclusion_scores: dict[str, float]
    selectivity: float
    metadata: dict = {}


class ExcludeManyResponse(BaseModel):
    gene_a: str
    excluded_genes: list[str]
    lineage: list[str]
    total_ranked: int
    lines: list[ExcludeManyLine]
    lineage_distribution: dict[str, int] = {}


class LineageRankedLine(BaseModel):
    model_id: str
    name: str | None = None
    core_score: float
    lineage: str
    rank_within_lineage: int
    rank_global: int
    metadata: dict = {}


class LineageRankedResponse(BaseModel):
    gene: str
    lineages_returned: int
    total_scoreable: int
    lineage_unassigned_count: int
    lines: list[LineageRankedLine]


class MultiLineageRankedLine(BaseModel):
    model_id: str
    name: str | None = None
    joint_score: float
    scores: dict[str, float | None]
    limiting_gene: str | None = None
    lineage: str
    rank_within_lineage: int
    rank_global: int
    metadata: dict = {}


class MultiLineageRankedResponse(BaseModel):
    genes: list[str]
    lineages_returned: int
    total_scoreable: int
    lineage_unassigned_count: int
    lines: list[MultiLineageRankedLine]


class ExcludeLineageRankedLine(BaseModel):
    model_id: str
    name: str | None = None
    score_a: float
    score_b: float
    selectivity: float
    lineage: str
    rank_within_lineage: int
    rank_global: int
    metadata: dict = {}


class ExcludeLineageRankedResponse(BaseModel):
    gene_a: str
    gene_b: str
    lineages_returned: int
    total_scoreable: int
    lineage_unassigned_count: int
    lines: list[ExcludeLineageRankedLine]


class ExcludeManyLineageRankedLine(BaseModel):
    model_id: str
    name: str | None = None
    score_a: float
    exclusion_scores: dict[str, float]
    selectivity: float
    lineage: str
    rank_within_lineage: int
    rank_global: int
    metadata: dict = {}


class ExcludeManyLineageRankedResponse(BaseModel):
    gene_a: str
    excluded_genes: list[str]
    lineages_returned: int
    total_scoreable: int
    lineage_unassigned_count: int
    lines: list[ExcludeManyLineageRankedLine]


class CellLineRef(BaseModel):
    model_id: str
    name: str | None = None


class DetailResponse(BaseModel):
    gene: str
    ensg: str
    cell_line: CellLineRef
    rank: int
    total: int
    score: float
    tier: str
    n_layers: int
    driver_alteration: bool
    p_mutation: float | None = None
    p_fusion: float | None = None
    has_cna_alteration: bool
    expression_level: float | None = None
    proteomics_level: float | None = None
    expression_by_source: dict[str, float | None] = {}
    proteomics_by_source: dict[str, float | None] = {}
    metadata: dict = {}
    rna_alternatives: list = []


class HealthResponse(BaseModel):
    status: str = "ok"
