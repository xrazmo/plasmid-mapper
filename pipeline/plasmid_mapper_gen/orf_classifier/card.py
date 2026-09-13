from .blast_search import best_hit_against_db, best_hit_against_db_diamond


def search_card(
    query_faa_path, db_prefix, molecule, min_identity, min_coverage,
    use_diamond=False, diamond_db_prefix=None, threads=4,
):
    """CARD is distributed as either a protein (translated homolog models)
    or nucleotide (gene sequences) FASTA depending on the release artifact
    downloaded -- molecule (from databases.detect_molecule()) picks the
    matching BLAST program: blastp for a protein database, tblastn (protein
    query vs six-frame-translated nucleotide database) for a nucleotide one.

    diamond has no tblastn equivalent, so it is only ever used for a
    protein-molecule database (use_diamond is only honored when
    molecule == "prot"); a nucleotide CARD release always uses tblastn
    regardless of whether diamond is installed.
    """
    if molecule == "prot" and use_diamond:
        return best_hit_against_db_diamond(
            query_faa_path, diamond_db_prefix, "CARD", min_identity, min_coverage, threads
        )
    program = "blastp" if molecule == "prot" else "tblastn"
    return best_hit_against_db(
        query_faa_path, db_prefix, "CARD", program, min_identity, min_coverage
    )
