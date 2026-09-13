# plasmid-mapper-gen

Generates `ref_data.js` and `pl_data.js` for [plasmid_mapper](../mapper.html)
from a query plasmid FASTA and a set of reference plasmid FASTAs — no more
hand-editing those files.

## Install

```bash
conda env create -f environment.yml
conda activate plasmid-mapper
```

This installs Prokka, NCBI BLAST+, `diamond`, and the `plasmid-mapper-gen`
command. `diamond` is used automatically whenever it is on `PATH` for any
protein-vs-protein database search (CARD/VFDB/BacMet/UniProt, whenever
that database's FASTA is protein) — it has no equivalent for a protein
query against a nucleotide subject, so ISfinder (always nucleotide) and
any nucleotide-molecule CARD/VFDB/BacMet release still use BLAST+'s
`tblastn` regardless. If `diamond` is missing, the pipeline falls back to
`blastp`/`tblastn` for everything with no functional difference in
results, just slower — this matters most for UniProt/SwissProt, by far
the largest of the five databases.

## Reference databases

You need local copies of up to five protein/nucleotide databases used to
classify each predicted ORF, checked in this priority order (a higher-
priority database's hit always wins, regardless of relative identity/
coverage): **CARD > ISfinder > VFDB > BacMet > UniProt/SwissProt**. Point
`--db-dir` at a directory containing any of the following (missing ones
are simply skipped — ORFs that would have matched them fall through to
the next database, or end up "hypothetical" if none match).

These files are large (UniProt alone is ~300MB decompressed) and are
never committed to this repository. The easiest way to populate a
`--db-dir` is the included download script, which fetches all five from
their real upstream sources and is safe to re-run later to refresh or
replace any of them:

```bash
python scripts/fetch_reference_databases.py --db-dir /path/to/dbs
# Re-fetch just one, e.g. after a new UniProt release:
python scripts/fetch_reference_databases.py --db-dir /path/to/dbs --only uniprot --force
```

| File in `--db-dir` | Database | Used for | Source (what the script fetches) |
|---|---|---|---|
| `card.fasta` | CARD | antibiotic resistance genes (`args`) | CARD v4.0.0 (pinned): `https://card.mcmaster.ca/download/0/broadstreet-v4.0.0.tar.bz2`, using its `protein_fasta_protein_homolog_model.fasta`. Pass `--card-version` to fetch a different release. |
| `isfinder.fasta` | ISfinder | insertion sequence elements (`isel`) | Nucleotide IS sequences from https://github.com/thanhleviet/ISfinder-sequences (`IS.fna`). |
| `vfdb.fasta` | VFDB (setA, curated core set) | virulence factors (`virulence`) | `https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz` |
| `bacmet.fasta` | BacMet | biocide/metal resistance genes (`biocidemetal`) | Experimentally-confirmed + predicted protein FASTA from http://bacmet.biomedicine.gu.se/, concatenated into one file. |
| `uniprot_sprot.fasta` | UniProt/SwissProt (reviewed) | general annotation fallback (`virulence`/`transposase`/`integrase`/`other`/`hypothetical`, by keyword match on the hit description — see note below) | Current release from `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`. |

UniProt's keyword-based "virulence" classification is a fallback, not the
primary signal — it only applies to an ORF that VFDB didn't classify
(missing `vfdb.fasta`, or no VFDB hit clearing `--min-identity`/
`--min-coverage`), since VFDB's setA is a curated core set rather than an
exhaustive one.

The first time you run the CLI against a raw FASTA, it builds a BLAST index
next to it automatically (`makeblastdb`); subsequent runs reuse that index.

## Usage

```bash
plasmid-mapper-gen \
  --query kpc33_plasmid.fasta --query-id KPC33_p1 \
  --reference ref1.fasta --reference ref2.fasta \
  --db-dir /path/to/dbs \
  --out-dir ./output
```

This writes `output/ref_data.js` (ORF annotations for the query plasmid)
and, if at least one `--reference` was given, `output/pl_data.js` (pairwise
BLASTN comparison rings against each reference). Copy both into
`plasmid_mapper/js/` to replace the existing files (or pass `--append` to
merge a new query plasmid into an existing `ref_data.js` in `--out-dir`
without discarding other plasmids already in it).

Each `--query`/`--reference` FASTA file must contain exactly one sequence.
Query IDs may only contain letters, digits, `.`, `_`, `-` (no `$`, which is
reserved as the query/subject separator in `pl_data.js`'s keys).

### Per-database identity/coverage thresholds

Each database has its own identity/coverage threshold an ORF's best hit
must clear to be classified via that tier (otherwise it falls through to
the next one). Defaults live in
[`plasmid_mapper_gen/db_thresholds.yaml`](plasmid_mapper_gen/db_thresholds.yaml):

| Database | min_identity | min_coverage |
|---|---|---|
| CARD | 70 | 70 |
| ISfinder | 70 | 70 |
| VFDB | 70 | 70 |
| BacMet | 70 | 70 |
| UniProt/SwissProt | 40 | 50 |

CARD/ISfinder/VFDB/BacMet stay strict, since a false-positive resistance/
virulence/IS/biocide-metal call is worse than a missed one. UniProt is
deliberately looser — it's the last fallback before "hypothetical", so a
weaker but real homology hit there is more useful than no annotation.

Override any subset of these with `--db-thresholds path/to/file.yaml`
(same format as the bundled file — you only need to include the
databases you want to change; anything omitted keeps its default):

```yaml
uniprot: { min_identity: 30, min_coverage: 40 }
```

### Useful flags

- `--db-thresholds PATH`: override per-database thresholds (see above).
- `--blastn-task megablast|blastn|dc-megablast` (default `megablast`): use
  plain `blastn` or `dc-megablast` instead of `megablast` if your reference
  plasmids are more distantly related than same-species comparisons.
- `--threads N`

## Running the tests

```bash
pip install -e ".[test]"
pytest tests/
```
