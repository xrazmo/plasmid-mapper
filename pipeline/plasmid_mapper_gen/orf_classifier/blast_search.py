import csv
import io
from dataclasses import dataclass

from ..external_tools import run

_OUTFMT_COLUMNS = [
    "sseqid", "stitle", "pident", "length", "qstart", "qend", "qlen",
]
_OUTFMT = "6 " + " ".join(_OUTFMT_COLUMNS)


@dataclass
class DbHit:
    dbname: str
    subject_id: str
    description: str
    pident: float
    coverage: float  # percent of query ORF covered by this HSP


def _parse_best_hit(stdout: str, dbname: str, min_identity: float, min_coverage: float):
    """Shared row-parsing logic for both BLAST+ and diamond output, which
    both produce the same tab-separated -outfmt/--outfmt 6 columns
    (_OUTFMT_COLUMNS) and both sort rows by bitscore descending within a
    query -- so in either case the first row is already the best hit,
    and there is at most one query here so no need to group by qseqid.
    """
    reader = csv.DictReader(
        io.StringIO(stdout), fieldnames=_OUTFMT_COLUMNS, delimiter="\t"
    )
    for row in reader:
        qlen = int(row["qlen"])
        coverage = 100.0 * (int(row["qend"]) - int(row["qstart"]) + 1) / qlen
        pident = float(row["pident"])
        if pident >= min_identity and coverage >= min_coverage:
            return DbHit(
                dbname=dbname,
                subject_id=row["sseqid"],
                description=row["stitle"],
                pident=pident,
                coverage=coverage,
            )
        return None  # first row is the best hit; if it fails thresholds, none will pass
    return None


def best_hit_against_db(
    query_fasta_path: str,
    db_prefix: str,
    dbname: str,
    program: str,
    min_identity: float,
    min_coverage: float,
):
    """Run `program` (blastp/tblastn) of a single-ORF-protein FASTA against
    a BLAST database and return the best-scoring hit passing the identity/
    coverage thresholds, or None if nothing qualifies.
    """
    result = run(
        [
            program,
            "-query", query_fasta_path,
            "-db", db_prefix,
            "-outfmt", _OUTFMT,
            "-max_target_seqs", "1",
        ],
        error_context=f"{program} search against {dbname}",
    )
    return _parse_best_hit(result.stdout, dbname, min_identity, min_coverage)


def best_hit_against_db_diamond(
    query_fasta_path: str,
    dmnd_prefix: str,
    dbname: str,
    min_identity: float,
    min_coverage: float,
    threads: int = 4,
):
    """diamond blastp equivalent of best_hit_against_db(), used whenever a
    database's molecule is protein and diamond is available -- diamond has
    no tblastn equivalent, so nucleotide-molecule databases never call
    this (see databases.ensure_diamond_db()'s docstring).
    """
    result = run(
        [
            "diamond", "blastp",
            "-q", query_fasta_path,
            "-d", dmnd_prefix,
            "--outfmt", "6", *_OUTFMT_COLUMNS,
            "-k", "1",
            "--threads", str(threads),
            "--quiet",
        ],
        error_context=f"diamond blastp search against {dbname}",
    )
    return _parse_best_hit(result.stdout, dbname, min_identity, min_coverage)
