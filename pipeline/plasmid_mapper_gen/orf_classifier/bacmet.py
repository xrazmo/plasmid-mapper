from .blast_search import best_hit_against_db, best_hit_against_db_diamond


def search_bacmet(
    query_faa_path, db_prefix, molecule, min_identity, min_coverage,
    use_diamond=False, diamond_db_prefix=None, threads=4,
):
    """BacMet ships both predicted (protein) and experimentally confirmed
    (also protein, but a user could point --db-dir at a nucleotide export)
    FASTA variants -- molecule picks blastp vs tblastn accordingly, same
    reasoning as search_card().

    diamond (no tblastn equivalent) is only used when molecule == "prot".
    """
    if molecule == "prot" and use_diamond:
        return best_hit_against_db_diamond(
            query_faa_path, diamond_db_prefix, "biocide and metal resistance database",
            min_identity, min_coverage, threads,
        )
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path,
        db_prefix,
        "biocide and metal resistance database",
        program,
        min_identity,
        min_coverage,
    )
