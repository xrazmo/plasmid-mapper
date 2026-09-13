import re

_TOKEN_RE = re.compile(r"\d+|[A-Za-z*-]{2}")


def parse_btop_to_line_annot(btop: str, qstart: int, qend: int, sstrand: str):
    """Walk a BLAST BTOP string and return a list of {"v": <query coord>,
    "t": <mismatch type>} entries, one per mismatch/gap position, in the
    schema pl_data.js's MAP_DATA[...].ranges[...].line_annot expects.

    BTOP alternates two token kinds:
      - a run of matching bases, encoded as a plain integer
      - a single mismatch/gap column, encoded as a 2-character pair
        "<query_base><subject_base>" (using "-" on whichever side has a gap)

    Query coordinates in BLASTN's tabular output (qstart/qend) are always
    reported low-to-high regardless of alignment strand, but the BTOP
    token stream itself walks the alignment in the order it was computed,
    which corresponds to increasing query coordinate when the query is
    plus-strand relative to the subject, and *decreasing* query coordinate
    when the subject alignment is minus-strand (sstrand == "minus") --
    BLAST always reports qstart<qend but conceptually walks the query
    3'->5' in that case. We therefore start the cursor at qend and
    decrement for minus-strand HSPs, at qstart and increment otherwise.

    Type mapping (a judgment call documented here since the existing
    hand-curated data's "m"/"s"/"q" distinction renders identically today
    via Mismatch_COLOR and isn't independently specified anywhere):
      - substitution (both sides real bases, differ)   -> "m"
      - subject has a gap ("-") at this query position  -> "q"
      - query has a gap ("-") at this position           -> "s" (no query
        coordinate exists to annotate; cursor does not advance and no
        line_annot entry is emitted for this token)
    """
    tokens = _TOKEN_RE.findall(btop)
    forward = sstrand != "minus"
    qpos = qstart if forward else qend
    step = 1 if forward else -1

    annotations = []
    for token in tokens:
        if token.isdigit():
            qpos += step * int(token)
            continue

        query_base, subject_base = token[0], token[1]
        if query_base == "-":
            # Query has a gap here: no query coordinate to annotate, cursor
            # does not advance (this column consumes no query base).
            continue
        if subject_base == "-":
            annotations.append({"v": qpos, "t": "q"})
        else:
            annotations.append({"v": qpos, "t": "m"})
        qpos += step

    return annotations
