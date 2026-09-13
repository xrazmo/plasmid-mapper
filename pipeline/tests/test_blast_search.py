import subprocess

from plasmid_mapper_gen.orf_classifier import blast_search

# A representative outfmt-6 row shaped exactly like _OUTFMT_COLUMNS
# (sseqid stitle pident length qstart qend qlen), tab-separated. Both
# BLAST+ (-outfmt 6 ...) and diamond (--outfmt 6 ...) produce this same
# column shape and both sort by bitscore descending, so parsing must
# treat their stdout identically.
_ROW = "VFG000676\thlyA hemolysin\t92.5\t100\t1\t100\t100\n"


def _fake_completed_process(stdout):
    return subprocess.CompletedProcess(args=["x"], returncode=0, stdout=stdout, stderr="")


def test_best_hit_against_db_parses_blast_row(monkeypatch):
    monkeypatch.setattr(
        blast_search, "run", lambda cmd, error_context: _fake_completed_process(_ROW)
    )
    hit = blast_search.best_hit_against_db(
        "query.faa", "db_prefix", "VFDB", "blastp", min_identity=70.0, min_coverage=70.0
    )
    assert hit is not None
    assert hit.subject_id == "VFG000676"
    assert hit.description == "hlyA hemolysin"
    assert hit.pident == 92.5
    assert hit.coverage == 100.0


def test_best_hit_against_db_diamond_parses_same_row_shape(monkeypatch):
    # The key guarantee this test locks in: diamond's outfmt-6 output
    # parses identically to BLAST+'s, since best_hit_against_db_diamond
    # shares _parse_best_hit with best_hit_against_db rather than
    # duplicating the coverage/threshold logic.
    monkeypatch.setattr(
        blast_search, "run", lambda cmd, error_context: _fake_completed_process(_ROW)
    )
    hit = blast_search.best_hit_against_db_diamond(
        "query.faa", "dmnd_prefix", "VFDB", min_identity=70.0, min_coverage=70.0, threads=4
    )
    assert hit is not None
    assert hit.subject_id == "VFG000676"
    assert hit.description == "hlyA hemolysin"
    assert hit.pident == 92.5
    assert hit.coverage == 100.0


def test_best_hit_against_db_diamond_invokes_diamond_blastp(monkeypatch):
    captured = {}

    def fake_run(cmd, error_context):
        captured["cmd"] = cmd
        return _fake_completed_process(_ROW)

    monkeypatch.setattr(blast_search, "run", fake_run)
    blast_search.best_hit_against_db_diamond(
        "query.faa", "dmnd_prefix", "VFDB", min_identity=70.0, min_coverage=70.0, threads=8
    )
    assert captured["cmd"][:2] == ["diamond", "blastp"]
    assert "-q" in captured["cmd"] and "query.faa" in captured["cmd"]
    assert "-d" in captured["cmd"] and "dmnd_prefix" in captured["cmd"]
    assert "--threads" in captured["cmd"] and "8" in captured["cmd"]


def test_below_threshold_returns_none(monkeypatch):
    low_identity_row = "VFG000676\thlyA hemolysin\t10.0\t100\t1\t100\t100\n"
    monkeypatch.setattr(
        blast_search, "run", lambda cmd, error_context: _fake_completed_process(low_identity_row)
    )
    hit = blast_search.best_hit_against_db(
        "query.faa", "db_prefix", "VFDB", "blastp", min_identity=70.0, min_coverage=70.0
    )
    assert hit is None


def test_no_hits_returns_none(monkeypatch):
    monkeypatch.setattr(
        blast_search, "run", lambda cmd, error_context: _fake_completed_process("")
    )
    hit = blast_search.best_hit_against_db(
        "query.faa", "db_prefix", "VFDB", "blastp", min_identity=70.0, min_coverage=70.0
    )
    assert hit is None
