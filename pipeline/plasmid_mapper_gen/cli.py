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
from .pairwise.blastn_runner import run_pairwise_blastn
from .prokka_runner import run_prokka
from .assembler.contig_ref import build_contig_ref_entry
from .assembler.map_data import build_map_data_entries
from .assembler.js_writer import (
    merge_contig_ref,
    write_pl_data_js,
    write_ref_data_js,
)
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
            "Directory containing card.fasta, bacmet.fasta, isfinder.fasta, "
            "uniprot_sprot.fasta (any may be omitted)"
        ),
    )
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--min-identity", type=float, default=70.0)
    parser.add_argument("--min-coverage", type=float, default=70.0)
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
    return parser


def run_pipeline(config: RunConfig) -> None:
    validate_plasmid_id(config.query_id)
    os.makedirs(config.out_dir, exist_ok=True)

    available_optional = check_tools_on_path()
    if not available_optional.get("diamond"):
        print(
            "Note: diamond not found on PATH; using blastp/tblastn "
            "(slower, but no functional difference in results).",
            file=sys.stderr,
        )

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

    orf_records = []
    for i, orf in enumerate(orfs):
        query_faa_path = os.path.join(config.out_dir, f"_orf_{i}.faa")
        with open(query_faa_path, "w") as f:
            f.write(f">{orf.locus_tag}\n{orf.protein_seq}\n")

        card_hit = (
            search_card(query_faa_path, db_prefixes["card"], config.min_identity, config.min_coverage)
            if "card" in db_prefixes
            else None
        )
        bacmet_hit = (
            search_bacmet(query_faa_path, db_prefixes["bacmet"], config.min_identity, config.min_coverage)
            if "bacmet" in db_prefixes
            else None
        )
        isfinder_hit = (
            search_isfinder(query_faa_path, db_prefixes["isfinder"], config.min_identity, config.min_coverage)
            if "isfinder" in db_prefixes
            else None
        )
        uniprot_hit = (
            search_uniprot(query_faa_path, db_prefixes["uniprot"], config.min_identity, config.min_coverage)
            if "uniprot" in db_prefixes
            else None
        )
        os.remove(query_faa_path)

        classification = classify_orf(card_hit, bacmet_hit, isfinder_hit, uniprot_hit)
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
        )
        map_data = build_map_data_entries(config.query_id, comparisons)
        pl_data_path = os.path.join(config.out_dir, "pl_data.js")
        write_pl_data_js(map_data, pl_data_path)
        print(f"Wrote {pl_data_path}", file=sys.stderr)
    else:
        print(
            "No --reference FASTAs given; skipping pairwise BLASTN "
            "(pl_data.js not written).",
            file=sys.stderr,
        )


def main(argv=None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    config = RunConfig(
        query_fasta=args.query,
        query_id=args.query_id,
        reference_fastas=args.references,
        db_dir=args.db_dir,
        out_dir=args.out_dir,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        blastn_task=args.blastn_task,
        threads=args.threads,
        append=args.append,
    )
    try:
        run_pipeline(config)
    except PipelineError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
