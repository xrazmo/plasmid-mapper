from .blast_search import best_hit_against_db, best_hit_against_db_diamond

# Ordered by specificity: checked top-to-bottom, first match wins. A
# maintainable single list rather than scattered inline string checks, so
# adding/tuning keywords doesn't require touching classification logic.
#
# The "virulence" rule here is a fallback, not the primary virulence
# signal: merge.py's classify_orf() now checks a dedicated VFDB tier
# ahead of UniProt, so this only ever fires for an ORF that VFDB missed
# (no vfdb.fasta in --db-dir, or its best VFDB hit didn't clear
# --min-identity/--min-coverage). VFDB's setA is a deliberately curated
# core set, not exhaustive, so a same-ORF UniProt keyword match at that
# point is a genuine second-best signal rather than noise VFDB was meant
# to replace -- removing it would silently downgrade such ORFs to "other"
# any time vfdb.fasta happens to be missing, which is worse than today's
# behavior, not better.
_KEYWORD_RULES = [
    ("transposase", ["transposase", "transposon"]),
    ("integrase", ["integrase", "integron"]),
    (
        "virulence",
        [
            "virulence", "toxin", "adhesin", "hemolysin", "haemolysin",
            "invasin", "fimbrial", "siderophore",
        ],
    ),
]


def search_uniprot(
    query_faa_path, db_prefix, molecule, min_identity, min_coverage,
    use_diamond=False, diamond_db_prefix=None, threads=4,
):
    """UniProt/SwissProt is always protein. Since it is by far the largest
    of the 5 databases (500k+ sequences), it is the search this benefits
    from diamond the most -- but the dispatch still follows molecule, for
    consistency with the other search_* functions and in case a
    differently-prepared --db-dir ever points uniprot at something else.
    """
    if molecule == "prot" and use_diamond:
        return best_hit_against_db_diamond(
            query_faa_path, diamond_db_prefix, "UniProt/SwissProt",
            min_identity, min_coverage, threads,
        )
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path,
        db_prefix,
        "UniProt/SwissProt",
        program,
        min_identity,
        min_coverage,
    )


def classify_by_keyword(description: str) -> str:
    """Classify a UniProt hit's product description into one of
    transposase/integrase/virulence/other/hypothetical by keyword match.

    Falls back to "hypothetical" only when the description itself indicates
    an unannotated protein; a confident but otherwise-unclassified hit is
    "other" rather than "hypothetical" (a real, described protein that just
    doesn't fall into any of the tracked categories).
    """
    lowered = description.lower()
    if "hypothetical protein" in lowered or not description.strip():
        return "hypothetical"
    for orf_type, keywords in _KEYWORD_RULES:
        if any(kw in lowered for kw in keywords):
            return orf_type
    return "other"
