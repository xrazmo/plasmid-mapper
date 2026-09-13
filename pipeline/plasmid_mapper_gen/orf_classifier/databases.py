import os
from dataclasses import dataclass

from ..external_tools import run
from ..utils.errors import PipelineError

_NUCLEOTIDE_CHARS = set("ACGTUNRYSWKMBDHV-")


def detect_molecule(fasta_path: str) -> str:
    """Sniff whether a FASTA file is nucleotide or protein by inspecting
    its first sequence's character composition.

    Database identity does not reliably imply molecule type: CARD, for
    instance, is distributed in both nucleotide (gene sequences) and
    protein (translated homolog models) FASTA files depending on which
    release artifact a user downloads, and hardcoding an assumption here
    previously caused every CARD search to silently return zero hits
    (blastp against a nucleotide database matches nothing) -- discovered
    when a real KPC-33 plasmid's Prokka-annotated bla_KPC ORF, confirmed
    present via Prokka's own generic annotation, was misclassified as
    "hypothetical" because the CARD FASTA in use was nucleotide.
    """
    with open(fasta_path) as f:
        f.readline()  # header
        sample = ""
        while len(sample) < 200:
            line = f.readline()
            if not line or line.startswith(">"):
                break
            sample += line.strip()
    if not sample:
        raise PipelineError(f"Could not read a sequence from {fasta_path}")
    nucleotide_fraction = sum(1 for c in sample.upper() if c in _NUCLEOTIDE_CHARS) / len(sample)
    return "nucl" if nucleotide_fraction > 0.95 else "prot"


@dataclass
class ReferenceDatabase:
    name: str  # human-readable, used verbatim as `dbname` in output
    fasta_path: str
    molecule: str  # "prot" or "nucl" -- see detect_molecule()


def ensure_blast_db(db: ReferenceDatabase, index_dir: str) -> str:
    """Ensure a BLAST index for `db` exists under index_dir, building it
    from the raw FASTA on first use if needed. Returns the BLAST DB path
    (without extension) to pass to -db.
    """
    if not os.path.exists(db.fasta_path):
        raise PipelineError(
            f"Reference database FASTA not found for {db.name}: {db.fasta_path}\n"
            "See pipeline/README.md for how to obtain/prepare this database."
        )

    os.makedirs(index_dir, exist_ok=True)
    db_prefix = os.path.join(index_dir, os.path.basename(db.fasta_path))
    ext = ".phr" if db.molecule == "prot" else ".nhr"
    # Newer BLAST+ may build version-5 DBs with a different extension;
    # check both to avoid rebuilding an already-valid index unnecessarily.
    already_built = os.path.exists(db_prefix + ext) or os.path.exists(
        db_prefix + ".00" + ext
    )
    if not already_built:
        run(
            [
                "makeblastdb",
                "-in", db.fasta_path,
                "-dbtype", db.molecule,
                "-out", db_prefix,
                "-title", db.name,
            ],
            error_context=f"Building BLAST database for {db.name}",
        )
    return db_prefix


def ensure_diamond_db(db: ReferenceDatabase, index_dir: str, threads: int = 4) -> str:
    """Ensure a diamond index for `db` exists under index_dir, building it
    from the raw FASTA on first use if needed. Returns the diamond DB
    path (without the .dmnd extension, which `diamond` appends itself) to
    pass to -d/--db.

    Only ever called for molecule == "prot" databases: diamond has no
    tblastn equivalent (protein query vs six-frame-translated nucleotide
    subject), so a nucleotide-molecule database (ISfinder always; CARD/
    VFDB/BacMet sometimes, depending on release) must keep using BLAST+'s
    tblastn regardless of whether diamond is installed. Callers are
    responsible for checking db.molecule == "prot" before calling this.
    """
    if not os.path.exists(db.fasta_path):
        raise PipelineError(
            f"Reference database FASTA not found for {db.name}: {db.fasta_path}\n"
            "See pipeline/README.md for how to obtain/prepare this database."
        )

    os.makedirs(index_dir, exist_ok=True)
    db_prefix = os.path.join(index_dir, os.path.basename(db.fasta_path) + ".dmnd_idx")
    if not os.path.exists(db_prefix + ".dmnd"):
        run(
            [
                "diamond", "makedb",
                "--in", db.fasta_path,
                "--db", db_prefix,
                "--threads", str(threads),
                "--quiet",
            ],
            error_context=f"Building diamond database for {db.name}",
        )
    return db_prefix


def build_registry(db_dir: str) -> dict:
    """Build the standard 5-database registry from a directory convention:
    <db_dir>/card.fasta, isfinder.fasta, vfdb.fasta, bacmet.fasta,
    uniprot_sprot.fasta. pipeline/scripts/fetch_reference_databases.py
    populates a directory in exactly this shape.

    Any of the five may be omitted (missing file); classify_orf() treats an
    absent database as "no hit" for that tier rather than failing the run,
    so a user without e.g. an ISfinder FASTA can still get partial results.

    Molecule type (nucleotide vs protein) is auto-detected per file via
    detect_molecule() rather than assumed from the database's identity --
    see that function's docstring for why a fixed assumption is unsafe.
    """
    names = {
        "card": "CARD",
        "isfinder": "ISFinder",
        "vfdb": "virulence factor database (VFDB)",
        "bacmet": "biocide and metal resistance database",
        "uniprot": "UniProt/SwissProt",
    }
    registry = {}
    for key, filename in [
        ("card", "card.fasta"),
        ("isfinder", "isfinder.fasta"),
        ("vfdb", "vfdb.fasta"),
        ("bacmet", "bacmet.fasta"),
        ("uniprot", "uniprot_sprot.fasta"),
    ]:
        fasta_path = os.path.join(db_dir, filename)
        if os.path.exists(fasta_path):
            registry[key] = ReferenceDatabase(
                names[key], fasta_path, detect_molecule(fasta_path)
            )
    return registry
