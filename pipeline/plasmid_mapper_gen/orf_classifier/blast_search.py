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

    BLAST's default output ordering is by bitscore descending within a
    query, so the first row is already the best hit; there is at most one
    query here so no need to group by qseqid.
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
    reader = csv.DictReader(
        io.StringIO(result.stdout), fieldnames=_OUTFMT_COLUMNS, delimiter="\t"
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
