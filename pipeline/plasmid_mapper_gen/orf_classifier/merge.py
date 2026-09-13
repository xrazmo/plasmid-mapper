import re
from dataclasses import dataclass

from .uniprot import classify_by_keyword

# CARD's FASTA headers are not consistently formatted across release
# artifacts. Two real formats confirmed by inspection:
#   1. "ARO:3004623|ID:3322|Name:AAC(3)-IId|NCBI:EU022314.1" (seen in an
#      older/alternately-packaged CARD copy)
#   2. "gb|ABS70977.1|ARO:3004623|AAC(3)-IId [Escherichia coli]" (the
#      actual format of protein_fasta_protein_homolog_model.fasta from
#      CARD's own official v4.0.0 download archive, confirmed by
#      downloading and inspecting it directly -- this is what
#      scripts/fetch_reference_databases.py's users will actually have)
# Both patterns are tried, first match wins.
_CARD_NAME_RE = re.compile(r"Name:([^|]+)")
_CARD_ARO_GENE_RE = re.compile(r"ARO:\d+\|([^\[|]+?)\s*(?:\[|$)")

# ISfinder/UniProt headers have no structured short-name field to extract
# (unlike CARD's ARO-tagged format), so their descriptions are just capped
# at a fixed length instead -- a plain safety net so the ORF popover/label
# text stays readable rather than showing a full raw BLAST subject title.
_ISFINDER_DESCRIPTION_MAX_LEN = 40
_UNIPROT_DESCRIPTION_MAX_LEN = 40

# UniProt/SwissProt's stitle is a full BLAST-style FASTA header, e.g.
# "sp|P62590|INT2_ECOLX Integrase/recombinase OS=Escherichia coli OX=562
# GN=intA PE=3 SV=1" -- the "sp|ACCESSION|NAME_ORG " prefix and the
# trailing "OS=.../OX=.../GN=.../PE=.../SV=..." annotation fields aren't
# the actual product description, but they used to eat most of the
# character budget before truncation, leaving almost nothing of the real
# text (e.g. "sp|P62590|INT2_ECOLX Integrase/recombinase OS=Escherichia
# c…" showed barely more than the accession/organism). Stripping both
# ends first means the truncation budget is spent on the description
# itself.
_UNIPROT_PREFIX_RE = re.compile(r"^\w{2}\|[^|]+\|\S+\s+")
_UNIPROT_SUFFIX_RE = re.compile(r"\s+OS=.*$")


@dataclass
class OrfClassification:
    type: str
    dbname: str
    refprotien: str
    idty: float
    cov: float
    dscr: str


def _card_short_name(description: str) -> str:
    """CARD's FASTA header is a pipe-delimited identifier string, not a
    human-readable product description like BacMet/UniProt/ISfinder
    headers already are, and its exact shape varies by release artifact
    (see the two confirmed formats documented above _CARD_NAME_RE/
    _CARD_ARO_GENE_RE). Extract just the gene/allele name so labels in
    the rendered figure read "AAC(3)-IId" instead of the raw header.
    Falls back to the raw description if neither pattern matches (e.g. a
    differently-formatted CARD release), so this never produces an empty
    label.
    """
    match = _CARD_NAME_RE.search(description)
    if match:
        return match.group(1)
    match = _CARD_ARO_GENE_RE.search(description)
    if match:
        return match.group(1).strip()
    return description


def _truncate(description: str, max_len: int) -> str:
    """Cap a raw BLAST subject title at max_len characters so it stays
    readable in a popover/label, ending in an ellipsis if it was cut.
    Descriptions already at or under the cap are returned unchanged.
    """
    stripped = description.strip()
    if len(stripped) <= max_len:
        return stripped
    return stripped[: max_len - 1].rstrip() + "…"


def _isfinder_short_description(description: str) -> str:
    return _truncate(description, _ISFINDER_DESCRIPTION_MAX_LEN)


def _uniprot_short_description(description: str) -> str:
    stripped = _UNIPROT_PREFIX_RE.sub("", description)
    stripped = _UNIPROT_SUFFIX_RE.sub("", stripped)
    return _truncate(stripped, _UNIPROT_DESCRIPTION_MAX_LEN)


def classify_orf(card_hit=None, isfinder_hit=None, vfdb_hit=None, bacmet_hit=None, uniprot_hit=None):
    """Resolve a single ORF's category from up to 5 independent database
    hits via a strict priority cascade (not a cross-database bitscore
    comparison, which isn't meaningful across differently sized/composed
    databases).

    Priority, highest first: CARD (args) > ISfinder (isel) > VFDB
    (virulence) > BacMet (biocidemetal) > UniProt/SwissProt
    (keyword-classified) > none (hypothetical). A hit in a higher tier
    always wins over a lower tier regardless of relative identity/
    coverage, as long as it clears its own configured threshold —
    best_hit_against_db() already enforces the threshold, so any non-None
    hit passed in here has already qualified.

    Parameter order matches the cascade order itself (not e.g. the order
    databases were first added historically), so the signature stays
    self-documenting as tiers are added/reordered.
    """
    if card_hit is not None:
        return OrfClassification(
            type="args",
            dbname=card_hit.dbname,
            refprotien=card_hit.subject_id,
            idty=card_hit.pident,
            cov=card_hit.coverage,
            dscr=_card_short_name(card_hit.description),
        )
    if isfinder_hit is not None:
        return OrfClassification(
            type="isel",
            dbname=isfinder_hit.dbname,
            refprotien=isfinder_hit.subject_id,
            idty=isfinder_hit.pident,
            cov=isfinder_hit.coverage,
            dscr=_isfinder_short_description(isfinder_hit.description),
        )
    if vfdb_hit is not None:
        return OrfClassification(
            type="virulence",
            dbname=vfdb_hit.dbname,
            refprotien=vfdb_hit.subject_id,
            idty=vfdb_hit.pident,
            cov=vfdb_hit.coverage,
            dscr=vfdb_hit.description,
        )
    if bacmet_hit is not None:
        return OrfClassification(
            type="biocidemetal",
            dbname=bacmet_hit.dbname,
            refprotien=bacmet_hit.subject_id,
            idty=bacmet_hit.pident,
            cov=bacmet_hit.coverage,
            dscr=bacmet_hit.description,
        )
    if uniprot_hit is not None:
        return OrfClassification(
            type=classify_by_keyword(uniprot_hit.description),
            dbname=uniprot_hit.dbname,
            refprotien=uniprot_hit.subject_id,
            idty=uniprot_hit.pident,
            cov=uniprot_hit.coverage,
            dscr=_uniprot_short_description(uniprot_hit.description),
        )
    return OrfClassification(
        type="hypothetical",
        dbname="",
        refprotien="",
        idty=0.0,
        cov=0.0,
        dscr="hypothetical protein",
    )
