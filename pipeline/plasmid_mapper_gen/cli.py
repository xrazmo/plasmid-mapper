import argparse
import os
import sys

from .config import RunConfig
from .external_tools import check_tools_on_path
from .orf_classifier import databases as db_registry
from .orf_classifier.bacmet import search_bacmet
from .orf_classifier.card import search_card
from .orf_classifier.isfinder import search_isfinder
from .orf_classifier.merge import classify_orf
from .orf_classifier.uniprot import search_uniprot
from .orf_classifier.vfdb import search_vfdb
from .pairwise.blastn_runner import run_pairwise_blastn
from .prokka_runner import run_prokka
from .assembler.contig_ref import build_contig_ref_entry
from .assembler.map_data import build_map_data_entries
from .assembler.js_writer import (
    merge_contig_ref,
    write_pl_data_js,
    write_ref_data_js,
)
from .assembler.single_html import build_single_html
from .threshold_config import load_thresholds
from .utils.errors import PipelineError
from .utils.fasta import read_single_fasta_record, validate_plasmid_id


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plasmid-mapper-gen",
        description=(
            "Generate plasmid_mapper's ref_data.js/pl_data.js from a query "
            "plasmid FASTA and one or more reference plasmid FASTAs."
        ),
    )
    parser.add_argument("--query", required=True, help="Query plasmid FASTA file")
    parser.add_argument(
        "--query-id", required=True, help="Identifier for the query plasmid"
    )
    parser.add_argument(
        "--reference",
        action="append",
        default=[],
        dest="references",
        help="Reference plasmid FASTA file (repeatable)",
    )
    parser.add_argument(
        "--db-dir",
        required=True,
        help=(
            "Directory containing card.fasta, isfinder.fasta, vfdb.fasta, "
            "bacmet.fasta, uniprot_sprot.fasta (any may be omitted). "
            "See scripts/fetch_reference_databases.py to populate one."
        ),
    )
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument(
        "--db-thresholds",
        dest="db_thresholds_path",
        default=None,
        help=(
            "Path to a YAML file overriding per-database min_identity/"
            "min_coverage thresholds (see plasmid_mapper_gen/db_thresholds.yaml "
            "for the bundled defaults and format). May override just the "
            "databases you want to change; omitted databases keep the "
            "bundled default."
        ),
    )
    parser.add_argument(
        "--blastn-task",
        choices=["megablast", "blastn", "dc-megablast"],
        default="megablast",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--append",
        action="store_true",
        help="Merge into an existing ref_data.js in --out-dir instead of overwriting",
    )
    parser.add_argument(
        "--merge-gap-bp",
        type=int,
        default=10,
        help=(
            "Merge adjacent BLASTN HSPs separated by a query gap smaller than "
            "this (a small indel) into one comparison ring segment instead of "
            "rendering a visually noisy break between them (default: 10)"
        ),
    )
    parser.add_argument(
        "--single-html",
        action="store_true",
        help=(
            "Also write <out-dir>/plasmid_viewer.html: one self-contained "
            "HTML file with the project's own CSS/JS inlined (minified via "
            "esbuild, which must be on PATH -- `npm install -g esbuild`) "
            "and the generated plasmid data inlined. Open it directly, no "
            "server or other files needed. Third-party libraries (jQuery/"
            "Bootstrap/D3/etc.) still load from their CDNs."
        ),
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help=(
            "Directory containing this project's css/ and js/ (the same "
            "layout as the repo checkout with mapper.html). Required with "
            "--single-html; defaults to the plasmid_mapper_gen package's "
            "own repo checkout location if omitted."
        ),
    )
    return parser


def _default_project_root() -> str:
    """Default --project-root for --single-html: assumes this file is at
    <repo_root>/pipeline/plasmid_mapper_gen/cli.py, the layout of an
    editable checkout (matches how this project is actually installed --
    `pip install -e .` per environment.yml). A user running from a
    different layout (e.g. a non-editable install elsewhere) should pass
    --project-root explicitly.
    """
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_pipeline(config: RunConfig) -> None:
    validate_plasmid_id(config.query_id)
    os.makedirs(config.out_dir, exist_ok=True)

    available_optional = check_tools_on_path()
    use_diamond = bool(available_optional.get("diamond"))
    if use_diamond:
        print(
            "diamond found on PATH; using it for all protein-vs-protein "
            "database searches (CARD/VFDB/BacMet/UniProt when their "
            "molecule is protein). ISfinder, and any database whose "
            "molecule is nucleotide, still use BLAST+'s tblastn -- "
            "diamond has no equivalent for a protein query against a "
            "nucleotide subject.",
            file=sys.stderr,
        )
    else:
        print(
            "Note: diamond not found on PATH; using blastp/tblastn "
            "(slower, but no functional difference in results).",
            file=sys.stderr,
        )

    thresholds = load_thresholds(config.db_thresholds_path)

    query_record = read_single_fasta_record(config.query_fasta)
    qlen = len(query_record.seq)

    print(f"Running Prokka on {config.query_fasta} ...", file=sys.stderr)
    orfs = run_prokka(
        config.query_fasta, config.query_id, config.out_dir, config.threads
    )
    print(f"Prokka predicted {len(orfs)} ORFs.", file=sys.stderr)

    registry = db_registry.build_registry(config.db_dir)
    if not registry:
        print(
            f"Warning: no reference databases found in {config.db_dir}; "
            "all ORFs will be classified as hypothetical/unclassified.",
            file=sys.stderr,
        )
    index_dir = os.path.join(config.out_dir, "blastdb", "orf_classifier")
    db_prefixes = {
        key: db_registry.ensure_blast_db(db, index_dir)
        for key, db in registry.items()
    }
    diamond_prefixes = {
        key: db_registry.ensure_diamond_db(db, index_dir, config.threads)
        for key, db in registry.items()
        if use_diamond and db.molecule == "prot"
    }

    orf_records = []
    for i, orf in enumerate(orfs):
        query_faa_path = os.path.join(config.out_dir, f"_orf_{i}.faa")
        with open(query_faa_path, "w") as f:
            f.write(f">{orf.locus_tag}\n{orf.protein_seq}\n")

        card_min_identity, card_min_coverage = thresholds["card"]
        isfinder_min_identity, isfinder_min_coverage = thresholds["isfinder"]
        vfdb_min_identity, vfdb_min_coverage = thresholds["vfdb"]
        bacmet_min_identity, bacmet_min_coverage = thresholds["bacmet"]
        uniprot_min_identity, uniprot_min_coverage = thresholds["uniprot"]

        card_hit = (
            search_card(
                query_faa_path, db_prefixes["card"], registry["card"].molecule,
                card_min_identity, card_min_coverage,
                use_diamond=use_diamond, diamond_db_prefix=diamond_prefixes.get("card"),
                threads=config.threads,
            )
            if "card" in db_prefixes
            else None
        )
        isfinder_hit = (
            search_isfinder(query_faa_path, db_prefixes["isfinder"], isfinder_min_identity, isfinder_min_coverage)
            if "isfinder" in db_prefixes
            else None
        )
        vfdb_hit = (
            search_vfdb(
                query_faa_path, db_prefixes["vfdb"], registry["vfdb"].molecule,
                vfdb_min_identity, vfdb_min_coverage,
                use_diamond=use_diamond, diamond_db_prefix=diamond_prefixes.get("vfdb"),
                threads=config.threads,
            )
            if "vfdb" in db_prefixes
            else None
        )
        bacmet_hit = (
            search_bacmet(
                query_faa_path, db_prefixes["bacmet"], registry["bacmet"].molecule,
                bacmet_min_identity, bacmet_min_coverage,
                use_diamond=use_diamond, diamond_db_prefix=diamond_prefixes.get("bacmet"),
                threads=config.threads,
            )
            if "bacmet" in db_prefixes
            else None
        )
        uniprot_hit = (
            search_uniprot(
                query_faa_path, db_prefixes["uniprot"], registry["uniprot"].molecule,
                uniprot_min_identity, uniprot_min_coverage,
                use_diamond=use_diamond, diamond_db_prefix=diamond_prefixes.get("uniprot"),
                threads=config.threads,
            )
            if "uniprot" in db_prefixes
            else None
        )
        os.remove(query_faa_path)

        classification = classify_orf(
            card_hit=card_hit,
            isfinder_hit=isfinder_hit,
            vfdb_hit=vfdb_hit,
            bacmet_hit=bacmet_hit,
            uniprot_hit=uniprot_hit,
        )
        orf_records.append(
            {
                "id": i,
                "sidx": orf.sidx,
                "eidx": orf.eidx,
                "strand": orf.strand,
                "type": classification.type,
                "dbname": classification.dbname,
                "refprotien": classification.refprotien,
                "idty": classification.idty,
                "cov": classification.cov,
                "gap": 0,
                "mismatch": 0,
                "dscr": classification.dscr,
                "linkout": None,
            }
        )

    contig_ref_entry = build_contig_ref_entry(
        config.query_id,
        qlen,
        query_record.description or config.query_id,
        "",
        orf_records,
    )
    new_contig_ref = {config.query_id: contig_ref_entry}

    ref_data_path = os.path.join(config.out_dir, "ref_data.js")
    if config.append:
        new_contig_ref = merge_contig_ref(ref_data_path, new_contig_ref)
    write_ref_data_js(new_contig_ref, ref_data_path)
    print(f"Wrote {ref_data_path}", file=sys.stderr)

    if config.reference_fastas:
        print(
            f"Running pairwise BLASTN against {len(config.reference_fastas)} "
            "reference plasmid(s) ...",
            file=sys.stderr,
        )
        comparisons = run_pairwise_blastn(
            config.query_fasta,
            config.query_id,
            qlen,
            config.reference_fastas,
            config.out_dir,
            config.blastn_task,
            config.threads,
            config.merge_gap_bp,
        )
        map_data = build_map_data_entries(config.query_id, comparisons)
        pl_data_path = os.path.join(config.out_dir, "pl_data.js")
        write_pl_data_js(map_data, pl_data_path)
        print(f"Wrote {pl_data_path}", file=sys.stderr)
    else:
        map_data = {}
        print(
            "No --reference FASTAs given; skipping pairwise BLASTN "
            "(pl_data.js not written).",
            file=sys.stderr,
        )

    if config.single_html:
        project_root = config.project_root or _default_project_root()
        single_html_path = os.path.join(config.out_dir, "plasmid_viewer.html")
        print(f"Building self-contained {single_html_path} ...", file=sys.stderr)
        build_single_html(
            new_contig_ref,
            map_data,
            title=config.query_id,
            out_path=single_html_path,
            project_root=project_root,
        )
        print(f"Wrote {single_html_path}", file=sys.stderr)


def main(argv=None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    config = RunConfig(
        query_fasta=args.query,
        query_id=args.query_id,
        reference_fastas=args.references,
        db_dir=args.db_dir,
        out_dir=args.out_dir,
        db_thresholds_path=args.db_thresholds_path,
        blastn_task=args.blastn_task,
        threads=args.threads,
        append=args.append,
        merge_gap_bp=args.merge_gap_bp,
        single_html=args.single_html,
        project_root=args.project_root,
    )
    try:
        run_pipeline(config)
    except PipelineError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
