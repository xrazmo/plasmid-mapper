from .blast_search import best_hit_against_db

# Ordered by specificity: checked top-to-bottom, first match wins. A
# maintainable single list rather than scattered inline string checks, so
# adding/tuning keywords doesn't require touching classification logic.
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


def search_uniprot(query_faa_path, db_prefix, min_identity, min_coverage):
    return best_hit_against_db(
        query_faa_path,
        db_prefix,
        "UniProt/SwissProt",
        "blastp",
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
