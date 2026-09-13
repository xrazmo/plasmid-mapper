from plasmid_mapper_gen.orf_classifier.blast_search import DbHit
from plasmid_mapper_gen.orf_classifier.merge import classify_orf

CARD_HIT = DbHit("CARD", "42734", "MCR-9", 100.0, 100.0)
BACMET_HIT = DbHit(
    "biocide and metal resistance database", "BAC123", "mercury resistance protein", 95.0, 90.0
)
ISFINDER_HIT = DbHit("ISFinder", "ISKpn26", "IS6 family transposase", 99.0, 98.0)
UNIPROT_TRANSPOSASE = DbHit(
    "UniProt/SwissProt", "P12345", "IS3 family transposase", 80.0, 85.0
)
UNIPROT_INTEGRASE = DbHit("UniProt/SwissProt", "P23456", "integron integrase IntI1", 85.0, 90.0)
UNIPROT_VIRULENCE = DbHit("UniProt/SwissProt", "P34567", "hemolysin A", 90.0, 95.0)
UNIPROT_OTHER = DbHit("UniProt/SwissProt", "P45678", "WbuC", 100.0, 99.0)
UNIPROT_HYPOTHETICAL = DbHit(
    "UniProt/SwissProt", "P56789", "hypothetical protein", 40.0, 50.0
)


def test_card_wins_over_everything():
    result = classify_orf(
        card_hit=CARD_HIT,
        bacmet_hit=BACMET_HIT,
        isfinder_hit=ISFINDER_HIT,
        uniprot_hit=UNIPROT_TRANSPOSASE,
    )
    assert result.type == "args"
    assert result.dbname == "CARD"
    assert result.refprotien == "42734"


def test_card_beats_weak_uniprot_specifically():
    # The exact scenario named in the plan: a strong CARD hit must win even
    # against an otherwise-plausible-looking UniProt hit.
    result = classify_orf(card_hit=CARD_HIT, uniprot_hit=UNIPROT_TRANSPOSASE)
    assert result.type == "args"


def test_bacmet_wins_over_isfinder_and_uniprot():
    result = classify_orf(
        bacmet_hit=BACMET_HIT, isfinder_hit=ISFINDER_HIT, uniprot_hit=UNIPROT_VIRULENCE
    )
    assert result.type == "biocidemetal"
    assert result.dbname == "biocide and metal resistance database"


def test_isfinder_wins_over_uniprot():
    result = classify_orf(isfinder_hit=ISFINDER_HIT, uniprot_hit=UNIPROT_OTHER)
    assert result.type == "isel"
    assert result.dbname == "ISFinder"


def test_uniprot_transposase_keyword():
    result = classify_orf(uniprot_hit=UNIPROT_TRANSPOSASE)
    assert result.type == "transposase"


def test_uniprot_integrase_keyword():
    result = classify_orf(uniprot_hit=UNIPROT_INTEGRASE)
    assert result.type == "integrase"


def test_uniprot_virulence_keyword():
    result = classify_orf(uniprot_hit=UNIPROT_VIRULENCE)
    assert result.type == "virulence"


def test_uniprot_other_when_no_keyword_matches():
    result = classify_orf(uniprot_hit=UNIPROT_OTHER)
    assert result.type == "other"


def test_uniprot_hypothetical_description_classified_hypothetical():
    result = classify_orf(uniprot_hit=UNIPROT_HYPOTHETICAL)
    assert result.type == "hypothetical"


def test_no_hit_anywhere_is_hypothetical():
    result = classify_orf()
    assert result.type == "hypothetical"
    assert result.dscr == "hypothetical protein"
    assert result.dbname == ""
