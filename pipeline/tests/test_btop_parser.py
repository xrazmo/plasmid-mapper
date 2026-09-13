from plasmid_mapper_gen.pairwise.btop_parser import parse_btop_to_line_annot


def test_pure_match_no_annotations():
    # "50" = 50 matching bases, no mismatches at all.
    result = parse_btop_to_line_annot("50", qstart=100, qend=149, sstrand="plus")
    assert result == []


def test_single_substitution_plus_strand():
    # 10 matches, then a mismatch (query=A, subject=T), then 5 matches.
    # qstart=1 (BLAST 1-based) -> first mismatch sits at query position
    # qstart + 10 = 11.
    result = parse_btop_to_line_annot("10AT5", qstart=1, qend=16, sstrand="plus")
    assert result == [{"v": 11, "t": "m"}]


def test_multiple_substitutions_plus_strand():
    result = parse_btop_to_line_annot("12AT34GC5", qstart=100, qend=151, sstrand="plus")
    # first mismatch at 100+12=112, then 34 matches -> second mismatch at 112+1+34=147
    assert result == [{"v": 112, "t": "m"}, {"v": 147, "t": "m"}]


def test_subject_gap_is_type_q():
    # subject has "-" at this column: insertion in query relative to subject.
    result = parse_btop_to_line_annot("5A-3", qstart=10, qend=17, sstrand="plus")
    assert result == [{"v": 15, "t": "q"}]


def test_query_gap_is_skipped_no_query_coordinate():
    # query has "-": no query coordinate exists at a pure subject-insertion
    # column, so nothing is emitted and the cursor does not advance.
    result = parse_btop_to_line_annot("5-A3AT2", qstart=10, qend=19, sstrand="plus")
    # after "5": qpos=15; "-A" consumes no query base, qpos stays 15;
    # "3": qpos=18; "AT": mismatch at 18, qpos=19
    assert result == [{"v": 18, "t": "m"}]


def test_minus_strand_cursor_counts_down_from_qend():
    # Same token layout as test_single_substitution_plus_strand, but on the
    # minus strand: cursor starts at qend and decrements, so the mismatch
    # lands near the qend side instead of the qstart side.
    result = parse_btop_to_line_annot("10AT5", qstart=1, qend=16, sstrand="minus")
    # start at qend=16, minus 10 -> 6, mismatch emitted at 6, then -1 -> 5
    assert result == [{"v": 6, "t": "m"}]


def test_minus_strand_multiple_substitutions_ordering():
    result = parse_btop_to_line_annot(
        "12AT34GC5", qstart=100, qend=151, sstrand="minus"
    )
    # cursor starts at 151, -12 -> 139 (mismatch), -1 -> 138, -34 -> 104 (mismatch)
    assert result == [{"v": 139, "t": "m"}, {"v": 104, "t": "m"}]
