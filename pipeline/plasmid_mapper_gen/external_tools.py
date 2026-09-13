import shutil
import subprocess

from .utils.errors import PipelineError

REQUIRED_TOOLS = ["prokka", "blastn", "makeblastdb", "blastp", "tblastn"]
OPTIONAL_TOOLS = ["diamond", "blastx"]


def check_tools_on_path(required=REQUIRED_TOOLS, optional=OPTIONAL_TOOLS) -> dict:
    """Verify required external tools are on PATH before doing any work.

    Fails fast with a plain-English, actionable message rather than letting
    a missing-tool error surface deep inside a subprocess call after
    minutes of Prokka/BLAST work have already run.
    """
    missing_required = [t for t in required if shutil.which(t) is None]
    if missing_required:
        raise PipelineError(
            "Missing required tool(s): "
            + ", ".join(missing_required)
            + ".\nInstall them via the provided conda environment:\n"
            "    conda env create -f environment.yml\n"
            "    conda activate plasmid-mapper\n"
            "then re-run this command."
        )
    available_optional = {t: shutil.which(t) is not None for t in optional}
    return available_optional


def run(cmd: list, error_context: str) -> subprocess.CompletedProcess:
    """Run a subprocess command, raising PipelineError with the command's
    stderr on failure instead of a raw CalledProcessError traceback.
    """
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
    except FileNotFoundError as exc:
        raise PipelineError(
            f"{error_context}: command not found: {cmd[0]}"
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise PipelineError(
            f"{error_context} failed (exit code {exc.returncode}):\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr:\n{exc.stderr}"
        ) from exc
    return result
