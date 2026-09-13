def build_map_data_entries(query_id, subject_comparisons):
    """Assemble MAP_DATA["<qry>$<subject>"] entries from a list of
    SubjectComparison objects (see pairwise/ranges_builder.py).

    Key format "qry$subject" matches the existing hand-curated convention
    in pl_data.js exactly (e.g. "s257ECL_2$NZ_KY863418.1"), since query
    IDs are validated (utils/fasta.validate_plasmid_id) to never contain
    "$" themselves.
    """
    entries = {}
    for comparison in subject_comparisons:
        key = f"{query_id}${comparison.sseqid}"
        entries[key] = {
            "qseqid": query_id,
            "qlen": comparison.qlen,
            "sseqid": comparison.sseqid,
            "stitle": comparison.stitle,
            "ranges": comparison.ranges,
        }
    return entries
