from dataclasses import dataclass, field


@dataclass
class RunConfig:
    query_fasta: str
    query_id: str
    reference_fastas: list
    db_dir: str
    out_dir: str
    min_identity: float = 70.0
    min_coverage: float = 70.0
    blastn_task: str = "megablast"
    threads: int = 4
    append: bool = False
    merge_gap_bp: int = 10
