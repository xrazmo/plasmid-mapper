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


def parse_blastn_output(blastn_stdout: str, qlen: int, subject_titles: dict, merge_gap_bp: int = 10):
    """Parse blastn -outfmt "6 ...btop sstrand" stdout into one
    SubjectComparison per distinct sseqid, each carrying a `ranges` list
    matching pl_data.js's MAP_DATA[...] schema.

    subject_titles maps sseqid -> a human-readable title (e.g. from the
    subject FASTA's description line), used to populate `stitle` the same
    way the existing hand-curated MAP_DATA entries do
    (e.g. "NZ_KY863418.1 Enterobacter asburiae ... plasmid pOXA436 ...").

    Adjacent HSPs separated by a small gap (below merge_gap_bp) are merged
    into one ranges entry via merge_adjacent_hsps() below, since otherwise
    a 1-2bp indel renders as a visually noisy hairline break between two
    ring segments for what is effectively one continuous alignment.
    """
    reader = csv.DictReader(
        io.StringIO(blastn_stdout), fieldnames=_OUTFMT_COLUMNS, delimiter="\t"
    )
    comparisons: dict = {}
    raw_hsps: dict = {}
    for row in reader:
        sseqid = row["sseqid"]
        if sseqid not in comparisons:
            comparisons[sseqid] = SubjectComparison(
                qseqid=row["qseqid"],
                qlen=qlen,
                sseqid=sseqid,
                stitle=subject_titles.get(sseqid, sseqid),
            )
            raw_hsps[sseqid] = []

        qstart, qend = int(row["qstart"]), int(row["qend"])
        line_annot = parse_btop_to_line_annot(
            row["btop"], qstart, qend, row["sstrand"]
        )
        raw_hsps[sseqid].append(
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

    for sseqid, hsps in raw_hsps.items():
        comparisons[sseqid].ranges = merge_adjacent_hsps(hsps, merge_gap_bp)

    return list(comparisons.values())


def merge_adjacent_hsps(hsps: list, merge_gap_bp: int) -> list:
    """Merge consecutive HSPs (sorted by qstart) into one ranges entry when
    they represent a single, effectively-continuous alignment interrupted
    only by a small indel, and mark genuine unaligned stretches (a real
    gap no HSP covers) more clearly with a "gap"-type line_annot entry at
    each boundary.

    Two HSPs merge only if ALL of:
      - same strand,
      - the query gap (next.qstart - current.qend) is positive and below
        merge_gap_bp -- a negative gap means the HSPs overlap in query
        coordinates (e.g. a duplicated/repeated region aligning twice);
        that is NOT a small indel and must not be merged, since doing so
        would fabricate a self-overlapping span,
      - subject coordinates are also collinear and contiguous in the same
        direction as the query gap (plus strand: next.sstart - current.send
        is also positive and below merge_gap_bp; minus strand: the mirrored
        check). Without this, two HSPs close in query position but mapping
        to distant/reordered subject coordinates (a translocation, or two
        copies of a repeated element) would merge into one arc implying
        synteny that does not exist -- a correctness bug in the rendered
        figure, not just a cosmetic one.

    A gap at or above merge_gap_bp (or one that fails the strand/
    collinearity checks) is left as two separate ranges entries, with a
    "gap" line_annot marker inserted at each boundary coordinate so a
    genuine unaligned stretch reads differently from an ordinary mismatch
    tick in the rendered figure.
    """
    if not hsps:
        return []

    ordered = sorted(hsps, key=lambda h: h["qstart"])
    merged = [dict(ordered[0])]

    for hsp in ordered[1:]:
        current = merged[-1]
        query_gap = hsp["qstart"] - current["qend"]

        same_strand = hsp["sstrand"] == current["sstrand"]
        query_gap_mergeable = 0 < query_gap < merge_gap_bp

        if same_strand and hsp["sstrand"] != "minus":
            subject_gap = hsp["sstart"] - current["send"]
        elif same_strand:
            subject_gap = current["sstart"] - hsp["send"]
        else:
            subject_gap = None
        subject_gap_mergeable = subject_gap is not None and 0 < subject_gap < merge_gap_bp

        if same_strand and query_gap_mergeable and subject_gap_mergeable:
            current["qend"] = hsp["qend"]
            if hsp["sstrand"] != "minus":
                current["send"] = hsp["send"]
            else:
                current["sstart"] = hsp["sstart"]
            current["line_annot"] = current["line_annot"] + hsp["line_annot"]
            total_length = current["length"] + hsp["length"]
            current["pident"] = (
                current["pident"] * current["length"] + hsp["pident"] * hsp["length"]
            ) / total_length
            current["mismatch"] += hsp["mismatch"]
            current["length"] = total_length
            current["gapopen"] += hsp["gapopen"]
            current["bitscore"] += hsp["bitscore"]
        else:
            if query_gap >= merge_gap_bp:
                current["line_annot"].append({"v": current["qend"], "t": "gap"})
                hsp = dict(hsp)
                hsp["line_annot"] = [{"v": hsp["qstart"], "t": "gap"}] + hsp["line_annot"]
            merged.append(dict(hsp))

    return merged
