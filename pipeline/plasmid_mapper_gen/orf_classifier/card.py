from .blast_search import best_hit_against_db


def search_card(query_faa_path, db_prefix, min_identity, min_coverage):
    return best_hit_against_db(
        query_faa_path, db_prefix, "CARD", "blastp", min_identity, min_coverage
    )
