def build_contig_ref_entry(query_id, qlen, definition, organism, orf_records, annotations=None):
    """Assemble one Contig_ref[qryId] entry matching ref_data.js's schema.

    orf_records: list of dicts already shaped like ref_data.js's orfs[]
    entries (id, sidx, eidx, strand, type, dbname, refprotien, idty, cov,
    gap, mismatch, dscr, linkout) — built by the caller from a ProkkaOrf
    plus its OrfClassification.

    annotations: list of {sidx, eidx} "zoomed band" regions. Each is given
    a stable id (ann-0, ann-1, ...) so the frontend edits-overlay can
    reference/remove base annotations without relying on fragile
    coordinate-string matching. Defaults to a single band spanning the
    whole plasmid, matching the fallback pattern already used for very
    short reference contigs in the hand-curated data (e.g. "annotations:
    [{sidx: 0, eidx: qlen}]").
    """
    if annotations is None:
        annotations = [{"sidx": 0, "eidx": qlen}]
    annotations_with_ids = [
        {"id": f"ann-{i}", "sidx": ann["sidx"], "eidx": ann["eidx"]}
        for i, ann in enumerate(annotations)
    ]

    return {
        "accession": query_id,
        "qlen": qlen,
        "definition": definition,
        "organism": organism,
        "orfs": orf_records,
        "annotations": annotations_with_ids,
    }
