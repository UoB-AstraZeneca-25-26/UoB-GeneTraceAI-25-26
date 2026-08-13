"""
Shared identifier-ambiguity detectors. (C6 -- detector hygiene.)

Why this module exists
---------------------
The composite-identifier scan lived in Python (`re`) and its regression
assertion lived in SQL (`LIKE '%;%'`). They drifted, and the drift caused a
FALSE PASS: the single CCLE composite is `cnk3/ipcef1`, a slash, so the
assertion reported 0 composites while the scan reported 1. One delimiter set,
defined once, used by both languages -- and an equivalence test that fails the
build if they ever diverge again.

Delimiters, and why each is here:
    ';'    ProCan symbol row, UniProt multi-gene entries
    ','    some GEO annotation fields, HGNC alias/prev symbol lists
    '|'    older Ensembl / NCBI headers, HGNC gene_group
    '/'    HGNC alias fields, read-through loci (CNK3/IPCEF1), cell-line names
    ' + '  fusion notation
"""
from __future__ import annotations

import re

# The single source of truth. Everything else is derived from this.
DELIM_CHARS = ";,|/"
FUSION_DELIM = r"\s\+\s"

DELIMS = re.compile(f"[{re.escape(DELIM_CHARS)}]|{FUSION_DELIM}")
# DuckDB's regexp_matches uses RE2; the character class and ' + ' translate
# directly, but the escaping differs, so derive it rather than retype it.
DELIM_SQL_PATTERN = f"[{DELIM_CHARS}]| \\+ "


def is_composite(value) -> bool:
    """True if `value` names more than one entity."""
    return value is not None and bool(DELIMS.search(str(value)))


def n_entities(value) -> int:
    """How many entities `value` names (1 if not composite)."""
    return 1 if value is None else len(re.split(DELIMS, str(value)))


def composite_sql(col: str) -> str:
    """SQL predicate matching the SAME strings as `is_composite`."""
    return f"regexp_matches({col}, '{DELIM_SQL_PATTERN}')"


# --------------------------------------------------------------- equivalence
# Strings that must classify identically in Python and in SQL. The first two
# are the real cases; the rest are the delimiters and near-miss negatives.
FIXTURES = [
    ("h4c1;h4c2;h4c3", True),          # ProCan protein group
    ("cnk3/ipcef1", True),             # CCLE read-through -- the false pass
    ("braf", False),
    ("a,b", True),
    ("a|b", True),
    ("bcl2l2-pabpn1", False),          # hyphen is NOT a delimiter
    ("tp53 + egfr", True),             # fusion notation
    ("gene+other", False),             # '+' without spaces is not fusion
    ("", False),
]


def assert_detectors_agree(con) -> None:
    """Fail loudly if the Python and SQL detectors ever diverge.

    `con` is a DuckDB connection (any connection -- this touches no tables).
    """
    rows = ", ".join(f"('{s}')" for s, _ in FIXTURES if s != "")
    sql = con.execute(
        f"SELECT v, {composite_sql('v')} AS m FROM (VALUES {rows}) t(v)"
    ).df()
    sql_map = dict(zip(sql["v"], sql["m"].astype(bool)))
    bad = []
    for s, expected in FIXTURES:
        py = is_composite(s)
        if py != expected:
            bad.append(f"python: {s!r} -> {py}, expected {expected}")
        if s != "" and sql_map.get(s) != expected:
            bad.append(f"sql:    {s!r} -> {sql_map.get(s)}, expected {expected}")
    if bad:
        raise AssertionError(
            "composite-identifier detectors disagree with the fixtures:\n  "
            + "\n  ".join(bad))


# ------------------------------------------------------------ join guards
# GDSC drug names are NOT unique: 71 names carry more than one DRUG_ID
# (selumetinib -> 1062/1498/1736; JQ1 -> 1218/163/2172, spanning MGH and
# Sanger). Any join on name fans out silently, inflating row counts and
# double-weighting whichever compound was screened twice.
FORBIDDEN_JOIN_KEYS = {
    "gdsc_compounds": ("DRUG_NAME", "join on DRUG_ID; 71 names map to >1 id"),
    "cosmic_cna": ("cosmic_sample_name",
                   "98% of rows already resolve by name -- see C5"),
}


def assert_no_name_join(table: str, keys) -> None:
    """Raise if a join uses an identifier known to be non-unique."""
    if table in FORBIDDEN_JOIN_KEYS:
        col, why = FORBIDDEN_JOIN_KEYS[table]
        if col in set(keys):
            raise AssertionError(
                f"{table}.{col} is not a unique key -- {why}")


if __name__ == "__main__":
    import duckdb
    assert_detectors_agree(duckdb.connect())
    print("detector equivalence: PASS")
    for s, _ in FIXTURES:
        print(f"  {s!r:<24s} composite={is_composite(s)!s:<6s} "
              f"n_entities={n_entities(s)}")
