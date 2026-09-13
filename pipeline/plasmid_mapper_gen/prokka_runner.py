import os
from dataclasses import dataclass

from .external_tools import run
from .utils.errors import PipelineError


@dataclass
class ProkkaOrf:
    locus_tag: str
    sidx: int  # 0-based start, matching ref_data.js's existing convention
    eidx: int  # 0-based-exclusive end
    strand: int  # 1 or -1
    product: str
    protein_seq: str
    nucleotide_seq: str


def run_prokka(query_fasta: str, query_id: str, out_dir: str, threads: int = 4):
    """Run Prokka on the query plasmid FASTA and return its predicted ORFs.

    Prokka writes <prefix>.faa (protein), <prefix>.ffn (nucleotide gene
    sequences), and <prefix>.gff (coordinates/strand/product) into out_dir.
    We parse the .gff for coordinates/strand/product and cross-reference
    .faa/.ffn for sequences, keyed by locus tag.
    """
    prokka_dir = os.path.join(out_dir, "prokka")
    os.makedirs(prokka_dir, exist_ok=True)

    run(
        [
            "prokka",
            "--outdir", prokka_dir,
            "--prefix", query_id,
            "--locustag", query_id,
            "--cpus", str(threads),
            "--force",
            query_fasta,
        ],
        error_context=f"Prokka annotation of {query_fasta}",
    )

    gff_path = os.path.join(prokka_dir, f"{query_id}.gff")
    faa_path = os.path.join(prokka_dir, f"{query_id}.faa")
    ffn_path = os.path.join(prokka_dir, f"{query_id}.ffn")
    for path in (gff_path, faa_path, ffn_path):
        if not os.path.exists(path):
            raise PipelineError(
                f"Expected Prokka output not found: {path}\n"
                "Prokka may have failed silently; check the prokka/ log files "
                f"in {prokka_dir} for details."
            )

    protein_seqs = _read_fasta_by_id(faa_path)
    nucleotide_seqs = _read_fasta_by_id(ffn_path)

    orfs = []
    for locus_tag, sidx, eidx, strand, product in _parse_gff_features(gff_path):
        if locus_tag not in protein_seqs:
            # Non-CDS features (rRNA, tRNA, etc.) have no protein sequence;
            # skip them, they don't participate in ORF/database classification.
            continue
        orfs.append(
            ProkkaOrf(
                locus_tag=locus_tag,
                sidx=sidx,
                eidx=eidx,
                strand=strand,
                product=product,
                protein_seq=protein_seqs[locus_tag],
                nucleotide_seq=nucleotide_seqs.get(locus_tag, ""),
            )
        )
    return orfs


def _read_fasta_by_id(path: str) -> dict:
    from Bio import SeqIO

    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}


def _parse_gff_features(gff_path: str):
    """Yield (locus_tag, sidx, eidx, strand, product) for each CDS feature.

    GFF3 coordinates are 1-based inclusive; converted here to the 0-based
    convention already used throughout ref_data.js's orfs[].sidx/eidx.
    """
    with open(gff_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break  # Prokka appends the full sequence after this marker
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            start_1based, end_1based = int(cols[3]), int(cols[4])
            strand = 1 if cols[6] == "+" else -1
            attrs = dict(
                item.split("=", 1) for item in cols[8].split(";") if "=" in item
            )
            locus_tag = attrs.get("locus_tag") or attrs.get("ID", "")
            product = attrs.get("product", "hypothetical protein")
            yield locus_tag, start_1based - 1, end_1based, strand, product
