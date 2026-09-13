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


def require_esbuild():
    """Only called from the --single-html output path -- esbuild is not
    needed for ordinary ref_data.js/pl_data.js generation, so it's not in
    REQUIRED_TOOLS (which would force every user to install it even if
    they never use --single-html).
    """
    if shutil.which("esbuild") is None:
        raise PipelineError(
            "--single-html requires esbuild, which was not found on PATH.\n"
            "Install it (after activating the plasmid-mapper conda env):\n"
            "    npm install -g esbuild\n"
            "then re-run this command."
        )


def run(cmd: list, error_context: str) -> subprocess.CompletedProcess:
    """Run a subprocess command, raising PipelineError with the command's
    stderr on failure instead of a raw CalledProcessError traceback.

    Decodes stdout/stderr with errors="replace" rather than the default
    strict UTF-8 decoding -- confirmed necessary in practice: some
    reference databases (e.g. VFDB) contain non-UTF-8 bytes (a raw 0xA0
    non-breaking space byte, at least) in their FASTA headers, which
    makeblastdb echoes verbatim into its own warning/error output. With
    strict decoding, capturing that output crashes with a
    UnicodeDecodeError before this function ever gets a chance to inspect
    makeblastdb's actual exit code -- a third-party tool's messy input
    data corrupting a header comment shouldn't take down the whole
    pipeline via an unrelated encoding crash.
    """
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, errors="replace", check=True
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
