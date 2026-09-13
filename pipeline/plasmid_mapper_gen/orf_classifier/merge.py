import re
from dataclasses import dataclass

from .uniprot import classify_by_keyword

_CARD_NAME_RE = re.compile(r"Name:([^|]+)")


@dataclass
class OrfClassification:
    type: str
    dbname: str
    refprotien: str
    idty: float
    cov: float
    dscr: str


def _card_short_name(description: str) -> str:
    """CARD's FASTA header is a pipe-delimited ARO identifier string, e.g.
    "ARO:3004623|ID:3322|Name:AAC(3)-IId|NCBI:EU022314.1" -- not a
    human-readable product description like BacMet/UniProt/ISfinder
    headers already are. Extract just the gene/allele name so labels in
    the rendered figure read "AAC(3)-IId" instead of the raw header.
    Falls back to the raw description if the pattern doesn't match (e.g.
    a differently-formatted CARD release), so this never produces an
    empty label.
    """
    match = _CARD_NAME_RE.search(description)
    return match.group(1) if match else description


def classify_orf(card_hit=None, bacmet_hit=None, isfinder_hit=None, uniprot_hit=None):
    """Resolve a single ORF's category from up to 4 independent database
    hits via a strict priority cascade (not a cross-database bitscore
    comparison, which isn't meaningful across differently sized/composed
    databases).

    Priority, highest first: CARD (args) > BacMet (biocidemetal) >
    ISfinder (isel) > UniProt/SwissProt (keyword-classified) > none
    (hypothetical). A hit in a higher tier always wins over a lower tier
    regardless of relative identity/coverage, as long as it clears its own
    configured threshold — best_hit_against_db() already enforces the
    threshold, so any non-None hit passed in here has already qualified.
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
    if bacmet_hit is not None:
        return OrfClassification(
            type="biocidemetal",
            dbname=bacmet_hit.dbname,
            refprotien=bacmet_hit.subject_id,
            idty=bacmet_hit.pident,
            cov=bacmet_hit.coverage,
            dscr=bacmet_hit.description,
        )
    if isfinder_hit is not None:
        return OrfClassification(
            type="isel",
            dbname=isfinder_hit.dbname,
            refprotien=isfinder_hit.subject_id,
            idty=isfinder_hit.pident,
            cov=isfinder_hit.coverage,
            dscr=isfinder_hit.description,
        )
    if uniprot_hit is not None:
        return OrfClassification(
            type=classify_by_keyword(uniprot_hit.description),
            dbname=uniprot_hit.dbname,
            refprotien=uniprot_hit.subject_id,
            idty=uniprot_hit.pident,
            cov=uniprot_hit.coverage,
            dscr=uniprot_hit.description,
        )
    return OrfClassification(
        type="hypothetical",
        dbname="",
        refprotien="",
        idty=0.0,
        cov=0.0,
        dscr="hypothetical protein",
    )
