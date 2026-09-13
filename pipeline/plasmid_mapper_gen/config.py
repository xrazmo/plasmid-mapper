from dataclasses import dataclass, field


@dataclass
class RunConfig:
    query_fasta: str
    query_id: str
    reference_fastas: list
    db_dir: str
    out_dir: str
    db_thresholds_path: str = None
    blastn_task: str = "megablast"
    threads: int = 4
    append: bool = False
    merge_gap_bp: int = 10
    single_html: bool = False
    project_root: str = None
