import re

from .errors import PipelineError

_SAFE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def validate_plasmid_id(plasmid_id: str) -> None:
    """Reject plasmid IDs that would break the "qry$subject" MAP_DATA key
    convention or produce an invalid JS object key when emitted unquoted.

    "$" is the query/subject separator used throughout MAP_DATA keys
    (e.g. "s082Km_2$NZ_KY863418.1"), so a plasmid ID containing "$" would
    make those keys ambiguous to split. Quote characters would break the
    emitted `var X = {...}` JS literal.
    """
    if not plasmid_id:
        raise PipelineError("Plasmid ID must not be empty.")
    if "$" in plasmid_id:
        raise PipelineError(
            f"Plasmid ID '{plasmid_id}' must not contain '$' "
            "(reserved as the query/subject separator in MAP_DATA keys)."
        )
    if not _SAFE_ID_RE.match(plasmid_id):
        raise PipelineError(
            f"Plasmid ID '{plasmid_id}' must contain only letters, digits, "
            "'.', '_' or '-'."
        )


def read_single_fasta_record(fasta_path: str):
    """Read a FASTA file expected to contain exactly one record (a single
    plasmid sequence) and return it as a Bio.SeqRecord.

    Raises PipelineError with a clear message if the file is missing, empty,
    or contains more than one record (ambiguous which sequence is "the"
    plasmid).
    """
    from Bio import SeqIO

    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    except FileNotFoundError as exc:
        raise PipelineError(f"FASTA file not found: {fasta_path}") from exc

    if not records:
        raise PipelineError(f"No sequences found in FASTA file: {fasta_path}")
    if len(records) > 1:
        raise PipelineError(
            f"{fasta_path} contains {len(records)} sequences; expected exactly "
            "one plasmid sequence per file. Split multi-FASTA files first."
        )
    return records[0]
