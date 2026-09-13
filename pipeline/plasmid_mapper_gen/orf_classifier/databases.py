import os
from dataclasses import dataclass

from ..external_tools import run
from ..utils.errors import PipelineError


@dataclass
class ReferenceDatabase:
    name: str  # human-readable, used verbatim as `dbname` in output
    fasta_path: str
    molecule: str  # "prot" or "nucl"


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


def build_registry(db_dir: str) -> dict:
    """Build the standard 4-database registry from a directory convention:
    <db_dir>/card.fasta, bacmet.fasta, isfinder.fasta, uniprot_sprot.fasta.

    Any of the four may be omitted (missing file); classify_orf() treats an
    absent database as "no hit" for that tier rather than failing the run,
    so a user without e.g. an ISfinder FASTA can still get partial results.
    """
    candidates = {
        "card": ReferenceDatabase("CARD", os.path.join(db_dir, "card.fasta"), "prot"),
        "bacmet": ReferenceDatabase(
            "biocide and metal resistance database",
            os.path.join(db_dir, "bacmet.fasta"),
            "prot",
        ),
        "isfinder": ReferenceDatabase(
            "ISFinder", os.path.join(db_dir, "isfinder.fasta"), "nucl"
        ),
        "uniprot": ReferenceDatabase(
            "UniProt/SwissProt",
            os.path.join(db_dir, "uniprot_sprot.fasta"),
            "prot",
        ),
    }
    return {key: db for key, db in candidates.items() if os.path.exists(db.fasta_path)}
