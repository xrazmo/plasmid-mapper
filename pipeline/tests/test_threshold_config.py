import textwrap

from plasmid_mapper_gen.threshold_config import load_thresholds


def test_bundled_defaults_load():
    thresholds = load_thresholds()
    assert thresholds["card"] == (70.0, 70.0)
    assert thresholds["isfinder"] == (70.0, 70.0)
    assert thresholds["vfdb"] == (70.0, 70.0)
    assert thresholds["bacmet"] == (70.0, 70.0)
    # UniProt is deliberately looser than the other 4 -- the whole point
    # of this change ("reduce the threshold... more protein can be
    # annotated") -- so this is a real assertion, not a placeholder.
    assert thresholds["uniprot"] == (40.0, 50.0)


def test_partial_override_merges_over_defaults(tmp_path):
    override_path = tmp_path / "thresholds.yaml"
    override_path.write_text(
        textwrap.dedent(
            """
            uniprot: { min_identity: 30, min_coverage: 40 }
            """
        )
    )
    thresholds = load_thresholds(str(override_path))
    assert thresholds["uniprot"] == (30.0, 40.0)
    # Everything not mentioned in the override file keeps the bundled
    # default -- a user shouldn't have to repeat all 5 databases just to
    # tune one.
    assert thresholds["card"] == (70.0, 70.0)
    assert thresholds["isfinder"] == (70.0, 70.0)
    assert thresholds["vfdb"] == (70.0, 70.0)
    assert thresholds["bacmet"] == (70.0, 70.0)


def test_unknown_key_in_override_is_ignored(tmp_path):
    override_path = tmp_path / "thresholds.yaml"
    override_path.write_text(
        textwrap.dedent(
            """
            some_future_db: { min_identity: 10, min_coverage: 10 }
            card: { min_identity: 90, min_coverage: 90 }
            """
        )
    )
    thresholds = load_thresholds(str(override_path))
    assert thresholds["card"] == (90.0, 90.0)
    assert "some_future_db" not in thresholds
