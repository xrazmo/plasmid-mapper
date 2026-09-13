from .blast_search import best_hit_against_db


def search_bacmet(query_faa_path, db_prefix, molecule, min_identity, min_coverage):
    """BacMet ships both predicted (protein) and experimentally confirmed
    (also protein, but a user could point --db-dir at a nucleotide export)
    FASTA variants -- molecule picks blastp vs tblastn accordingly, same
    reasoning as search_card().
    """
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path,
        db_prefix,
        "biocide and metal resistance database",
        program,
        min_identity,
        min_coverage,
    )
