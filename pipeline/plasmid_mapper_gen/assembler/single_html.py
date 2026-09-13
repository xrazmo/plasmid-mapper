import importlib.resources
import os
import tempfile

from ..external_tools import require_esbuild, run
from .js_writer import build_var_assignment_js

# Project's own JS, concatenated in the same order mapper.html's <script>
# tags already use (minus mustache.js/templates.js, which are dead weight
# on that page -- confirmed unused by pl_mapper.js/pl_editor.js/
# pl_style_modal.js, loaded only for an unrelated page's popover feature).
# Nothing here relies on load ORDER at true top-level parse time (every
# cross-file reference happens inside $(document).ready(...)/
# $(function(){...}), which only fires after all scripts have loaded), but
# this order is kept anyway as the safe, already-tested default.
_PROJECT_JS_FILES = [
    "inline_style.js",
    "FileSaver.js",
    "pl_editor.js",
    "pl_mapper.js",
    "pl_style_modal.js",
]
_PROJECT_CSS_FILES = ["style.css", "sidebar.css"]


def _read_template() -> str:
    ref = importlib.resources.files("plasmid_mapper_gen.assembler") / "single_html_template.html"
    with importlib.resources.as_file(ref) as path:
        with open(path) as f:
            return f.read()


def _minify_js(source: str) -> str:
    """Minify JS via esbuild (a single global-scope script, no ES module
    syntax involved -- semantically identical to today's separate
    <script> tags, just concatenated first). Requires esbuild on PATH;
    call require_esbuild() before this to fail with an actionable message
    instead of a raw FileNotFoundError from run().
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        src_path = os.path.join(tmp_dir, "bundle.js")
        out_path = os.path.join(tmp_dir, "bundle.min.js")
        with open(src_path, "w") as f:
            f.write(source)
        run(
            ["esbuild", src_path, "--minify", f"--outfile={out_path}"],
            error_context="Minifying bundled JS with esbuild",
        )
        with open(out_path) as f:
            return f.read()


def render_single_html(
    title: str,
    inline_css: str,
    inline_js: str,
    contig_ref: dict,
    map_data: dict,
) -> str:
    """Substitute the single-HTML template's placeholders. Pure string
    templating, no esbuild/file-IO -- factored out from
    build_single_html() so it's unit-testable with a fake pre-minified JS
    string, independent of whether esbuild is installed.

    Uses plain str.replace() on distinct {{...}} tokens rather than
    str.format()/f-strings: the substituted CSS/JS content is full of
    literal '{' and '}' characters (JS blocks, CSS rules), which would
    collide with str.format()'s own brace syntax.
    """
    template = _read_template()
    data_js = (
        build_var_assignment_js("Contig_ref", contig_ref)
        + "\n"
        + build_var_assignment_js("MAP_DATA", map_data)
    )
    replacements = {
        "{title}": title,
        "{inline_css}": inline_css,
        "{inline_js}": inline_js,
        "{data_js}": data_js,
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template


def build_single_html(
    contig_ref: dict,
    map_data: dict,
    title: str,
    out_path: str,
    project_root: str,
) -> None:
    """Assemble one self-contained HTML file: the project's own CSS/JS
    inlined (JS minified via esbuild), the generated plasmid data
    inlined, and CDN <script>/<link> tags for jQuery/Bootstrap/Popper/D3/
    Font Awesome left as-is (per explicit decision: only the project's
    own assets get inlined, not third-party libraries -- keeps the build
    simple and these libraries need no minification work from us).

    project_root must contain css/ and js/ (the same layout as the repo
    checkout containing mapper.html) -- passed explicitly rather than
    assumed from cli.py's own file location, since plasmid-mapper-gen is
    pip-installable and could run from a layout where that assumption
    doesn't hold.
    """
    require_esbuild()

    css_dir = os.path.join(project_root, "css")
    js_dir = os.path.join(project_root, "js")

    css_parts = []
    for filename in _PROJECT_CSS_FILES:
        with open(os.path.join(css_dir, filename)) as f:
            css_parts.append(f.read())
    inline_css = "\n".join(css_parts)

    js_parts = []
    for filename in _PROJECT_JS_FILES:
        with open(os.path.join(js_dir, filename)) as f:
            js_parts.append(f.read())
    combined_js = "\n;\n".join(js_parts)
    inline_js = _minify_js(combined_js)

    html = render_single_html(title, inline_css, inline_js, contig_ref, map_data)
    with open(out_path, "w") as f:
        f.write(html)
