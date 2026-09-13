from .blast_search import best_hit_against_db


def search_vfdb(query_faa_path, db_prefix, molecule, min_identity, min_coverage):
    """VFDB setA (the curated core virulence-factor set) is distributed as
    a protein FASTA, but molecule (from databases.detect_molecule()) still
    picks blastp vs tblastn rather than assuming protein -- same reasoning
    as search_card(): a user pointing --db-dir at an unexpected molecule
    type should get a real search, not a silent zero-hit blastp-against-
    nucleotide mismatch.
    """
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path, db_prefix, "virulence factor database (VFDB)", program, min_identity, min_coverage
    )
