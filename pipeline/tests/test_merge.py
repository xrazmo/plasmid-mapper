from plasmid_mapper_gen.orf_classifier.blast_search import DbHit
from plasmid_mapper_gen.orf_classifier.merge import classify_orf

CARD_HIT = DbHit("CARD", "42734", "MCR-9", 100.0, 100.0)
CARD_HIT_RAW_ARO_HEADER = DbHit(
    "CARD",
    "3322",
    "ARO:3004623|ID:3322|Name:AAC(3)-IId|NCBI:EU022314.1",
    100.0,
    100.0,
)
CARD_HIT_OFFICIAL_DOWNLOAD_HEADER = DbHit(
    "CARD",
    "AUD37084.1",
    "gb|AUD37084.1|ARO:3005365|KPC-33 [Klebsiella pneumoniae] ",
    100.0,
    100.0,
)
BACMET_HIT = DbHit(
    "biocide and metal resistance database", "BAC123", "mercury resistance protein", 95.0, 90.0
)
ISFINDER_HIT = DbHit("ISFinder", "ISKpn26", "IS6 family transposase", 99.0, 98.0)
VFDB_HIT = DbHit("virulence factor database (VFDB)", "VFG000676", "hlyA hemolysin", 92.0, 96.0)
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
        isfinder_hit=ISFINDER_HIT,
        vfdb_hit=VFDB_HIT,
        bacmet_hit=BACMET_HIT,
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


def test_card_beats_vfdb_specifically():
    result = classify_orf(card_hit=CARD_HIT, vfdb_hit=VFDB_HIT)
    assert result.type == "args"


def test_isfinder_wins_over_bacmet_and_uniprot():
    # ISfinder now outranks BacMet (previous order had this reversed) --
    # this test replaces the old test_bacmet_wins_over_isfinder_and_uniprot,
    # whose assertion is false under the new priority order.
    result = classify_orf(
        bacmet_hit=BACMET_HIT, isfinder_hit=ISFINDER_HIT, uniprot_hit=UNIPROT_VIRULENCE
    )
    assert result.type == "isel"
    assert result.dbname == "ISFinder"


def test_isfinder_wins_over_vfdb_bacmet_and_uniprot():
    result = classify_orf(
        isfinder_hit=ISFINDER_HIT, vfdb_hit=VFDB_HIT, bacmet_hit=BACMET_HIT, uniprot_hit=UNIPROT_VIRULENCE
    )
    assert result.type == "isel"


def test_isfinder_wins_over_uniprot():
    result = classify_orf(isfinder_hit=ISFINDER_HIT, uniprot_hit=UNIPROT_OTHER)
    assert result.type == "isel"
    assert result.dbname == "ISFinder"


def test_vfdb_wins_over_bacmet_and_uniprot():
    # The key new-tier test: VFDB now sits between ISfinder and BacMet,
    # giving a real BLAST-based virulence call priority over both BacMet
    # and UniProt's weaker keyword-based virulence guess.
    result = classify_orf(vfdb_hit=VFDB_HIT, bacmet_hit=BACMET_HIT, uniprot_hit=UNIPROT_VIRULENCE)
    assert result.type == "virulence"
    assert result.dbname == "virulence factor database (VFDB)"
    assert result.refprotien == "VFG000676"


def test_vfdb_wins_over_uniprot_keyword_virulence():
    # Both VFDB and UniProt's keyword match would independently produce
    # type == "virulence" here -- assert the real VFDB hit's fields are
    # actually used, not UniProt's, guarding against a future refactor
    # that short-circuits on type equality instead of tier precedence.
    result = classify_orf(vfdb_hit=VFDB_HIT, uniprot_hit=UNIPROT_VIRULENCE)
    assert result.type == "virulence"
    assert result.refprotien == "VFG000676"
    assert result.dscr == "hlyA hemolysin"


def test_bacmet_wins_over_uniprot():
    # BacMet dropped below ISfinder/VFDB in the new order but still
    # outranks UniProt -- not covered by any surviving test after
    # test_isfinder_wins_over_bacmet_and_uniprot replaced the old
    # bacmet-vs-isfinder assertion above.
    result = classify_orf(bacmet_hit=BACMET_HIT, uniprot_hit=UNIPROT_VIRULENCE)
    assert result.type == "biocidemetal"
    assert result.dbname == "biocide and metal resistance database"


def test_uniprot_transposase_keyword():
    result = classify_orf(uniprot_hit=UNIPROT_TRANSPOSASE)
    assert result.type == "transposase"


def test_uniprot_integrase_keyword():
    result = classify_orf(uniprot_hit=UNIPROT_INTEGRASE)
    assert result.type == "integrase"


def test_no_vfdb_falls_through_to_uniprot_virulence_keyword():
    # No vfdb_hit passed at all -- simulates a missing vfdb.fasta in
    # --db-dir or a VFDB search that didn't clear the identity/coverage
    # threshold. UniProt's keyword-based virulence guess is the intended
    # fallback in that case (see the comment above _KEYWORD_RULES in
    # uniprot.py), not a hidden bug.
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


def test_card_dscr_extracts_short_name_from_raw_aro_header():
    # Real CARD FASTA headers are pipe-delimited ARO identifier strings,
    # not human-readable descriptions like BacMet/UniProt/ISfinder already
    # provide -- dscr should show just the gene name, not the raw header.
    result = classify_orf(card_hit=CARD_HIT_RAW_ARO_HEADER)
    assert result.dscr == "AAC(3)-IId"


def test_card_dscr_extracts_short_name_from_official_download_header():
    # CARD's actual official v4.0.0 download archive
    # (protein_fasta_protein_homolog_model.fasta, as fetched by
    # scripts/fetch_reference_databases.py) uses a DIFFERENT header format
    # than CARD_HIT_RAW_ARO_HEADER above -- confirmed by downloading and
    # inspecting the real file directly, which is what surfaced this case:
    # "gb|<accession>|ARO:<id>|<GeneName> [<organism>]", no "Name:" field
    # at all. Both formats must be handled, since a user following
    # pipeline/README.md's instructions will get this format, not the
    # other one.
    result = classify_orf(card_hit=CARD_HIT_OFFICIAL_DOWNLOAD_HEADER)
    assert result.dscr == "KPC-33"


def test_card_dscr_falls_back_to_raw_description_if_unparseable():
    unparseable_hit = DbHit("CARD", "999", "some other format entirely", 100.0, 100.0)
    result = classify_orf(card_hit=unparseable_hit)
    assert result.dscr == "some other format entirely"
