import os

from Bio import SeqIO

from ..external_tools import run
from .ranges_builder import BLASTN_OUTFMT, parse_blastn_output


def run_pairwise_blastn(
    query_fasta: str,
    query_id: str,
    qlen: int,
    reference_fastas: list,
    out_dir: str,
    blastn_task: str = "megablast",
    threads: int = 4,
):
    """Run BLASTN of the query plasmid against each reference plasmid FASTA
    in turn, returning a list of SubjectComparison objects (one or more per
    reference, since a reference FASTA could in principle contain multiple
    subject sequences).
    """
    blastdb_dir = os.path.join(out_dir, "blastdb")
    os.makedirs(blastdb_dir, exist_ok=True)

    all_comparisons = []
    for reference_fasta in reference_fastas:
        subject_titles = {
            rec.id: rec.description for rec in SeqIO.parse(reference_fasta, "fasta")
        }

        db_prefix = os.path.join(
            blastdb_dir, os.path.basename(reference_fasta)
        )
        run(
            [
                "makeblastdb",
                "-in", reference_fasta,
                "-dbtype", "nucl",
                "-out", db_prefix,
            ],
            error_context=f"Building BLAST database for {reference_fasta}",
        )

        result = run(
            [
                "blastn",
                "-query", query_fasta,
                "-db", db_prefix,
                "-task", blastn_task,
                "-outfmt", BLASTN_OUTFMT,
                "-num_threads", str(threads),
            ],
            error_context=f"BLASTN of {query_id} against {reference_fasta}",
        )
        comparisons = parse_blastn_output(result.stdout, qlen, subject_titles)
        all_comparisons.extend(comparisons)

    return all_comparisons
