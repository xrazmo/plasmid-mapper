# Plasmid Mapper

A D3.js tool for visualizing and comparing circular plasmids as
publication-ready figures. It renders a query plasmid as a circular map
annotated with ORFs (resistance genes, insertion sequences, transposons,
virulence factors, etc.) and overlays BLAST comparison rings against a
set of related reference plasmids — similar in spirit to tools like BRIG
or Proksee's ring viewer.

A companion Python pipeline (`pipeline/`) generates the map data
directly from FASTA files — running Prokka, classifying every predicted
ORF against five reference databases (CARD, ISfinder, VFDB, BacMet,
UniProt/SwissProt), and running pairwise BLASTN against each reference
plasmid — so there's no manual data editing involved.

## Try it

[`examples/KPC33_p1/kpc33_viewer.html`](examples/KPC33_p1/kpc33_viewer.html)
is a full worked example (a KPC-33 plasmid compared against two
reference plasmids), generated end-to-end by the pipeline. Open it in a
browser to see the tool in action.

## Features

- Circular ORF map with a curved or boxed category legend, colored by
  type (ARGs, insertion sequences, transposons, virulence factors,
  biocide/metal resistance, integrons, hypothetical, other).
- BLAST comparison rings against any number of reference plasmids, each
  with its own results-table row, mismatch/gap tick marks, and an
  independently pickable ring color.
- Interactive editing: draggable legends, zoom bands (magnify a region
  of the query ring into a second, zoomed-in ring), freeform and
  per-ORF labels, all persisted in the browser and exportable/importable
  as a JSON edits file.
- A Style modal for customizing font family/size per text element and
  color per ORF category.
- Export the current figure as SVG, PNG, or JPEG at a chosen resolution.

## Contents

- `mapper.html` — the circular plasmid/BLAST-ring viewer (generic
  template — starts with no plasmid data loaded; see `pipeline/` or
  `examples/` for real data to view).
- `js/pl_mapper.js` — core D3 rendering logic (rings, ORF arrows,
  legends, axis ticks).
- `js/pl_editor.js` — interactive editing layer (zoom bands, labels,
  legend dragging, ORF tooltips, edits import/export, image export).
- `js/pl_style_modal.js` — the font/color Style modal.
- `js/inline_style.js` — default color palettes for ORF categories and
  BLAST hit rings.
- `js/ref_data.js`, `js/pl_data.js` — the query plasmid's ORF
  annotations (`Contig_ref`) and pairwise BLASTN comparison data
  (`MAP_DATA`) that `mapper.html` loads. Ship empty by default;
  populated by the pipeline (or by hand, following the same shape).
- `js/FileSaver.js` — vendored client-side file-saving library (used for
  image and edits-file export).
- `css/style.css`, `css/sidebar.css` — styling (figure/table styling and
  the collapsible-sidebar app shell, respectively).
- `pipeline/` — the Python data-generation pipeline (see
  [`pipeline/README.md`](pipeline/README.md)).
- `examples/` — worked examples with real generated data.
- `index.html`, `html/aln.html` — a separate, unrelated linear alignment
  browser (not covered here).

## Usage

Open `mapper.html` in a browser (or serve the directory with any static
file server). If `js/ref_data.js`/`js/pl_data.js` have data in them,
select a query plasmid from the dropdown and click **Update**; choose
which BLAST subject hits to display as rings from the results table;
use the sidebar to adjust figure geometry, zoom bands, labels, and
style; export the figure from the Export image panel.

To generate real data instead of hand-editing `js/ref_data.js`, use the
pipeline (see [`pipeline/README.md`](pipeline/README.md)):

```bash
plasmid-mapper-gen \
  --query my_plasmid.fasta --query-id MyPlasmid \
  --reference ref1.fasta --reference ref2.fasta \
  --db-dir /path/to/dbs \
  --out-dir ./output \
  --single-html
```

`--single-html` writes one self-contained `output/plasmid_viewer.html`
with the viewer's own CSS/JS inlined and minified and the generated data
inlined — open it directly, no server or other files needed. Without
that flag, the pipeline writes `output/ref_data.js`/`output/pl_data.js`
to copy into `js/` instead.

## License

MIT — see [LICENSE](LICENSE).
