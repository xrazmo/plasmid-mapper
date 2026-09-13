from .blast_search import best_hit_against_db


def search_card(query_faa_path, db_prefix, molecule, min_identity, min_coverage):
    """CARD is distributed as either a protein (translated homolog models)
    or nucleotide (gene sequences) FASTA depending on the release artifact
    downloaded -- molecule (from databases.detect_molecule()) picks the
    matching BLAST program: blastp for a protein database, tblastn (protein
    query vs six-frame-translated nucleotide database) for a nucleotide one.
    """
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path, db_prefix, "CARD", program, min_identity, min_coverage
    )
