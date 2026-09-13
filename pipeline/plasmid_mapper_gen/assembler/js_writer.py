import json


def write_ref_data_js(contig_ref: dict, out_path: str):
    """Emit a ref_data.js-compatible file: `var Contig_ref = {...};`.

    Loaded via a plain <script> tag by mapper.html (no fetch/JSON.parse in
    the frontend), so the output must be valid JS, not just valid JSON —
    json.dumps produces JS-compatible object/array/string literals for the
    plain-data shapes used here, so this is a straight serialization.

    MCR_LOC is deliberately not emitted: confirmed dead in pl_mapper.js
    (declared, never read), fully superseded by Contig_ref[...].annotations.
    """
    with open(out_path, "w") as f:
        f.write("var Contig_ref = ")
        json.dump(contig_ref, f, indent=2)
        f.write(";\n")


def write_pl_data_js(map_data: dict, out_path: str):
    """Emit a pl_data.js-compatible file: `var MAP_DATA = {...};`."""
    with open(out_path, "w") as f:
        f.write("var MAP_DATA = ")
        json.dump(map_data, f, indent=2)
        f.write(";\n")


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
