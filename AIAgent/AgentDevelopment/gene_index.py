"""Local gene lookup index — the primary source for ``gene_alias_lookup``.

Loads ``reference/gene_lookup.parquet`` once into in-memory dicts keyed by
ENSG ID and by upper-cased HGNC symbol (aliases and previous symbols folded
in as secondary keys). Sub-millisecond lookups, no network round trip.

This addresses plan.md R1: Ensembl measured at 0% success during review, so
the 19,213-gene local table (same fields the tool returns) becomes the
primary source; Ensembl/HGNC become enrichment for genes outside the panel.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class GeneRecord:
    ensg_id: str
    symbol: str
    full_name: str | None
    previous_symbols: list[str]
    synonyms: list[str]


@dataclass
class GeneIndex:
    by_ensg: dict[str, GeneRecord] = field(default_factory=dict)
    by_symbol: dict[str, GeneRecord] = field(default_factory=dict)

    def lookup(self, query: str) -> GeneRecord | None:
        q = query.strip().upper()
        if q.startswith("ENSG"):
            return self.by_ensg.get(q.split(".")[0])
        return self.by_symbol.get(q)

    def __len__(self) -> int:
        return len(self.by_ensg)


def _split_list(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and value != value):  # NaN
        return []
    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if not text:
        return []
    return [v.strip() for v in text.split("|") if v.strip()]


def load_gene_index(path: Path) -> GeneIndex:
    """Load the parquet table into an in-memory ``GeneIndex``.

    Returns an empty (not an error) index if the file is missing, so a
    checkout without the reference data still starts — every lookup just
    falls through to the Ensembl/HGNC network sources instead.
    """
    index = GeneIndex()
    if not Path(path).exists():
        return index

    df = pd.read_parquet(path)
    for row in df.itertuples(index=False):
        ensg_id = str(row.ensg_id).strip()
        symbol = str(row.hgnc_symbol).strip() if row.hgnc_symbol else ""
        prev = _split_list(getattr(row, "prev_symbols", None))
        aliases = _split_list(getattr(row, "alias_symbols", None))
        full_name = getattr(row, "approved_name", None)
        full_name = str(full_name).strip() if full_name and str(full_name) != "nan" else None

        record = GeneRecord(
            ensg_id=ensg_id,
            symbol=symbol,
            full_name=full_name,
            previous_symbols=prev,
            synonyms=aliases,
        )

        index.by_ensg[ensg_id.upper()] = record
        if symbol:
            index.by_symbol.setdefault(symbol.upper(), record)
        for alias in (*prev, *aliases):
            index.by_symbol.setdefault(alias.upper(), record)

    return index
