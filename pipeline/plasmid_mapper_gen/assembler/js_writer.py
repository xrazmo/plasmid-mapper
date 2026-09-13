import json


def build_var_assignment_js(var_name: str, data: dict) -> str:
    """Build a `var <var_name> = <data as JS-compatible JSON>;` string --
    the single source of truth for this format, shared by the standalone
    ref_data.js/pl_data.js writers below and by the single-HTML assembler
    (assembler/single_html.py), which inlines the same text into a
    <script> block instead of a separate file.

    json.dumps produces JS-compatible object/array/string literals for
    the plain-data shapes used here, so this is a straight serialization,
    not a real JS-generation step.
    """
    return f"var {var_name} = {json.dumps(data, indent=2)};\n"


def write_ref_data_js(contig_ref: dict, out_path: str):
    """Emit a ref_data.js-compatible file: `var Contig_ref = {...};`.

    Loaded via a plain <script> tag by mapper.html (no fetch/JSON.parse in
    the frontend), so the output must be valid JS, not just valid JSON.

    MCR_LOC is deliberately not emitted: confirmed dead in pl_mapper.js
    (declared, never read), fully superseded by Contig_ref[...].annotations.
    """
    with open(out_path, "w") as f:
        f.write(build_var_assignment_js("Contig_ref", contig_ref))


def write_pl_data_js(map_data: dict, out_path: str):
    """Emit a pl_data.js-compatible file: `var MAP_DATA = {...};`."""
    with open(out_path, "w") as f:
        f.write(build_var_assignment_js("MAP_DATA", map_data))


def merge_contig_ref(existing_path: str, new_entries: dict) -> dict:
    """Load an existing ref_data.js-style file (if present) and merge new
    Contig_ref entries into it, so regenerating one plasmid doesn't discard
    others already in the file. Existing entries for the same query ID are
    overwritten by the new ones (re-running the pipeline for a plasmid is
    expected to replace its data).

    This is a minimal text-based extraction (not a JS parser): it evals the
    file's `Contig_ref` object via a narrow, explicit substitution, which is
    safe here because these files are pipeline-generated JSON-shaped data,
    not arbitrary script.
    """
    import os
    import re

    if not os.path.exists(existing_path):
        return new_entries

    with open(existing_path) as f:
        content = f.read()
    match = re.search(r"var\s+Contig_ref\s*=\s*(\{.*\});?\s*$", content, re.DOTALL)
    if not match:
        return new_entries

    existing = json.loads(match.group(1))
    existing.update(new_entries)
    return existing
