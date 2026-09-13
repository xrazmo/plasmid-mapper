from plasmid_mapper_gen.pairwise.ranges_builder import merge_adjacent_hsps


def make_hsp(qstart, qend, sstart, send, sstrand="plus", pident=100.0, length=None, mismatch=0, gapopen=0, bitscore=100.0):
    length = length if length is not None else (qend - qstart + 1)
    return {
        "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
        "pident": pident, "evalue": 0.0, "bitscore": bitscore, "mismatch": mismatch,
        "length": length, "gapopen": gapopen, "sstrand": sstrand, "line_annot": [],
    }


def test_small_gap_merges_into_one_range():
    hsps = [
        make_hsp(1, 100, 1, 100),
        make_hsp(103, 200, 103, 200),  # 2bp query gap, collinear subject
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 1
    assert result[0]["qstart"] == 1
    assert result[0]["qend"] == 200
    assert result[0]["sstart"] == 1
    assert result[0]["send"] == 200


def test_gap_exactly_at_threshold_does_not_merge():
    # gap == merge_gap_bp is NOT < merge_gap_bp (strict inequality), so
    # must stay separate: boundary/off-by-one case.
    hsps = [
        make_hsp(1, 100, 1, 100),
        make_hsp(110, 200, 110, 200),  # gap = 110-100 = 10 == threshold(10)
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 2


def test_gap_one_below_threshold_merges():
    hsps = [
        make_hsp(1, 100, 1, 100),
        make_hsp(109, 200, 109, 200),  # gap = 9, merge_gap_bp=10 -> 9 < 10, merges
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 1


def test_gap_above_threshold_stays_separate_with_gap_markers():
    hsps = [
        make_hsp(1, 100, 1, 100),
        make_hsp(150, 200, 150, 200),  # 50bp gap
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 2
    # boundary gap markers inserted on both sides of the real gap
    assert {"v": 100, "t": "gap"} in result[0]["line_annot"]
    assert {"v": 150, "t": "gap"} in result[1]["line_annot"]


def test_overlapping_hsps_do_not_merge():
    # negative query gap: HSPs overlap in query coordinates (e.g. a
    # duplicated/repeated region aligning twice) -- must not be merged as
    # if it were a small positive gap.
    hsps = [
        make_hsp(1, 100, 1, 100),
        make_hsp(90, 200, 90, 200),  # qstart(90) < current qend(100) -> overlap
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 2


def test_opposite_strand_adjacent_hsps_do_not_merge():
    hsps = [
        make_hsp(1, 100, 1, 100, sstrand="plus"),
        make_hsp(103, 200, 103, 200, sstrand="minus"),
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 2


def test_non_collinear_subject_despite_small_query_gap_does_not_merge():
    # Query gap is small (2bp) but subject coordinates jump far away --
    # a translocation or repeated element, not a simple indel. Must not
    # merge, since doing so would fabricate a contiguous synteny that
    # doesn't exist.
    hsps = [
        make_hsp(1, 100, 1, 100),
        make_hsp(103, 200, 5000, 5097),  # subject gap = 5000-100 = 4900, way above threshold
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 2


def test_minus_strand_merge_uses_mirrored_subject_check():
    # On minus strand, subject coordinates decrease as query increases;
    # collinearity check must mirror this (current.sstart - hsp.send).
    hsps = [
        make_hsp(1, 100, 200, 300, sstrand="minus"),
        make_hsp(103, 200, 100, 197, sstrand="minus"),  # subject continues downward
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 1
    assert result[0]["sstart"] == 100
    assert result[0]["send"] == 300


def test_pident_is_length_weighted_average_on_merge():
    hsps = [
        make_hsp(1, 100, 1, 100, pident=100.0, length=100),
        make_hsp(103, 152, 103, 152, pident=90.0, length=50),
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 1
    expected = (100.0 * 100 + 90.0 * 50) / 150
    assert abs(result[0]["pident"] - expected) < 1e-9


def test_empty_input_returns_empty_list():
    assert merge_adjacent_hsps([], merge_gap_bp=10) == []


def test_unsorted_input_gets_sorted_by_qstart():
    hsps = [
        make_hsp(103, 200, 103, 200),
        make_hsp(1, 100, 1, 100),
    ]
    result = merge_adjacent_hsps(hsps, merge_gap_bp=10)
    assert len(result) == 1
    assert result[0]["qstart"] == 1
