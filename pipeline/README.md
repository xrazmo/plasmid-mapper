# plasmid-mapper-gen

Generates `ref_data.js` and `pl_data.js` for [plasmid_mapper](../mapper.html)
from a query plasmid FASTA and a set of reference plasmid FASTAs — no more
hand-editing those files.

## Install

```bash
conda env create -f environment.yml
conda activate plasmid-mapper
```

This installs Prokka, NCBI BLAST+, and the `plasmid-mapper-gen` command.
(`diamond` is not required — BLAST's own `blastp`/`tblastn` are used
instead; if you have `diamond` installed separately it is not currently
used, but nothing breaks either way.)

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

### Useful flags

- `--min-identity` / `--min-coverage` (default 70/70): thresholds an ORF's
  best database hit must clear to be classified via that database, instead
  of falling through to the next tier.
- `--blastn-task megablast|blastn|dc-megablast` (default `megablast`): use
  plain `blastn` or `dc-megablast` instead of `megablast` if your reference
  plasmids are more distantly related than same-species comparisons.
- `--threads N`

## Running the tests

```bash
pip install -e ".[test]"
pytest tests/
```
