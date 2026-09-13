from plasmid_mapper_gen.assembler.single_html import render_single_html

_SAMPLE_CONTIG_REF = {
    "TestPlasmid": {
        "accession": "TestPlasmid",
        "qlen": 100,
        "definition": "TestPlasmid",
        "organism": "",
        "orfs": [],
        "annotations": [],
    }
}


def test_render_single_html_substitutes_all_placeholders():
    html = render_single_html(
        title="TestPlasmid",
        inline_css=".foo { color: red; }",
        inline_js="function foo() { return 1; }",
        contig_ref=_SAMPLE_CONTIG_REF,
        map_data={},
    )
    # No leftover {token} placeholders from the template.
    for token in ("{title}", "{inline_css}", "{inline_js}", "{data_js}"):
        assert token not in html

    assert "TestPlasmid" in html
    assert ".foo { color: red; }" in html
    assert "function foo() { return 1; }" in html
    assert "var Contig_ref" in html
    assert "var MAP_DATA" in html


def test_render_single_html_has_no_local_script_or_link_refs():
    # A single-file output must not reference any local project path --
    # only CDN URLs and inline <style>/<script> blocks.
    html = render_single_html(
        title="TestPlasmid",
        inline_css="",
        inline_js="",
        contig_ref=_SAMPLE_CONTIG_REF,
        map_data={},
    )
    assert 'src="js/' not in html
    assert 'href="css/' not in html
    assert 'href="./css/' not in html
    assert "mustache.js" not in html
    assert "templates.js" not in html


def test_render_single_html_includes_footer():
    html = render_single_html(
        title="TestPlasmid",
        inline_css="",
        inline_js="",
        contig_ref=_SAMPLE_CONTIG_REF,
        map_data={},
    )
    assert "app-footer" in html
    assert "plasmid-mapper" in html.lower()


def test_render_single_html_inlines_map_data():
    html = render_single_html(
        title="TestPlasmid",
        inline_css="",
        inline_js="",
        contig_ref=_SAMPLE_CONTIG_REF,
        map_data={"TestPlasmid$Ref1": {"qseqid": "TestPlasmid", "sseqid": "Ref1"}},
    )
    assert "Ref1" in html
    assert "TestPlasmid$Ref1" in html
