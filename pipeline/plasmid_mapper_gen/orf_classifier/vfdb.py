from .blast_search import best_hit_against_db, best_hit_against_db_diamond


def search_vfdb(
    query_faa_path, db_prefix, molecule, min_identity, min_coverage,
    use_diamond=False, diamond_db_prefix=None, threads=4,
):
    """VFDB setA (the curated core virulence-factor set) is distributed as
    a protein FASTA, but molecule (from databases.detect_molecule()) still
    picks blastp vs tblastn rather than assuming protein -- same reasoning
    as search_card(): a user pointing --db-dir at an unexpected molecule
    type should get a real search, not a silent zero-hit blastp-against-
    nucleotide mismatch.

    diamond (no tblastn equivalent) is only used when molecule == "prot".
    """
    if molecule == "prot" and use_diamond:
        return best_hit_against_db_diamond(
            query_faa_path, diamond_db_prefix, "virulence factor database (VFDB)",
            min_identity, min_coverage, threads,
        )
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path, db_prefix, "virulence factor database (VFDB)", program, min_identity, min_coverage
    )
