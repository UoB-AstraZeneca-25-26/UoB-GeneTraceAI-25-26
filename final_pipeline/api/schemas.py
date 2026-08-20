"""Request / response models for the final_pipeline ranking API."""

from pydantic import BaseModel, Field


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


<<<<<<< Updated upstream
class DetailResponse(BaseModel):
    gene: str
    ensg_id: str
    model_id: str
    name: str | None = None
    score: float
    rank: int
    total_lines: int
    lineage: dict
=======
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
    metadata: dict = {}
    rna_alternatives: list = []
>>>>>>> Stashed changes


class HealthResponse(BaseModel):
    status: str = "ok"
