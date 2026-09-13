from .blast_search import best_hit_against_db


def search_isfinder(query_faa_path, db_prefix, min_identity, min_coverage):
    """ISfinder is a nucleotide (insertion sequence) database, while the ORF
    query is a protein sequence — tblastn (protein query vs a
    six-frame-translated nucleotide database) is the correct tool here, not
    blastx (which expects a nucleotide query).
    """
    return best_hit_against_db(
        query_faa_path, db_prefix, "ISFinder", "tblastn", min_identity, min_coverage
    )
