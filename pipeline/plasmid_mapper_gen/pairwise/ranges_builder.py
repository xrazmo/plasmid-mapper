import csv
import io
from dataclasses import dataclass, field

from .btop_parser import parse_btop_to_line_annot

_OUTFMT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "btop", "sstrand",
]
BLASTN_OUTFMT = "6 " + " ".join(_OUTFMT_COLUMNS)


@dataclass
class SubjectComparison:
    qseqid: str
    qlen: int
    sseqid: str
    stitle: str
    ranges: list = field(default_factory=list)


def parse_blastn_output(blastn_stdout: str, qlen: int, subject_titles: dict):
    """Parse blastn -outfmt "6 ...btop sstrand" stdout into one
    SubjectComparison per distinct sseqid, each carrying a `ranges` list
    (one entry per raw HSP) matching pl_data.js's MAP_DATA[...] schema.

    subject_titles maps sseqid -> a human-readable title (e.g. from the
    subject FASTA's description line), used to populate `stitle` the same
    way the existing hand-curated MAP_DATA entries do
    (e.g. "NZ_KY863418.1 Enterobacter asburiae ... plasmid pOXA436 ...").
    """
    reader = csv.DictReader(
        io.StringIO(blastn_stdout), fieldnames=_OUTFMT_COLUMNS, delimiter="\t"
    )
    comparisons: dict = {}
    for row in reader:
        sseqid = row["sseqid"]
        if sseqid not in comparisons:
            comparisons[sseqid] = SubjectComparison(
                qseqid=row["qseqid"],
                qlen=qlen,
                sseqid=sseqid,
                stitle=subject_titles.get(sseqid, sseqid),
            )

        qstart, qend = int(row["qstart"]), int(row["qend"])
        line_annot = parse_btop_to_line_annot(
            row["btop"], qstart, qend, row["sstrand"]
        )
        comparisons[sseqid].ranges.append(
            {
                "qstart": qstart,
                "qend": qend,
                "sstart": int(row["sstart"]),
                "send": int(row["send"]),
                "pident": float(row["pident"]),
                "evalue": float(row["evalue"]),
                "bitscore": float(row["bitscore"]),
                "mismatch": int(row["mismatch"]),
                "length": int(row["length"]),
                "gapopen": int(row["gapopen"]),
                "sstrand": row["sstrand"],
                "line_annot": line_annot,
            }
        )
    return list(comparisons.values())
