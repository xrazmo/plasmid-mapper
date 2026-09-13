#!/usr/bin/env python3
"""Fetches the 5 reference databases plasmid-mapper-gen classifies ORFs
against (CARD, ISfinder, VFDB, BacMet, UniProt/SwissProt) into a
--db-dir, matching the plain-filename convention
plasmid_mapper_gen.orf_classifier.databases.build_registry() expects
(card.fasta, isfinder.fasta, vfdb.fasta, bacmet.fasta,
uniprot_sprot.fasta).

Deliberately a standalone script, not a plasmid-mapper-gen CLI
subcommand: these files are large (UniProt alone is ~300MB decompressed)
and are never meant to be committed to the repo or bundled with the
installed package -- this script is how a user (re-)populates a --db-dir
on their own machine, on their own schedule.

Usage:
    python fetch_reference_databases.py --db-dir /path/to/dbs
    python fetch_reference_databases.py --db-dir /path/to/dbs --force
    python fetch_reference_databases.py --db-dir /path/to/dbs --only uniprot,vfdb
    python fetch_reference_databases.py --db-dir /path/to/dbs --card-version 4.1.0

Notes:
- BacMet's two source files are served over plain HTTP (not HTTPS) --
  this is a property of BacMet's own server, not a mistake here. Some
  networks/proxies block or warn on plain HTTP; if BacMet fails to fetch,
  that's the most likely cause.
- CARD is pinned to a specific version (default 4.0.0) rather than
  chasing "latest" the way UniProt's URL does, since CARD has no
  "current release" URL and its tarball's internal file layout is not
  guaranteed stable across versions -- re-fetching a different version is
  a deliberate, explicit choice via --card-version, and this script fails
  loudly (rather than silently picking a wrong file) if the expected
  protein_fasta_protein_homolog_model.fasta member isn't found after
  extraction.
"""

import argparse
import gzip
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request

CARD_URL_TEMPLATE = "https://card.mcmaster.ca/download/0/broadstreet-v{version}.tar.bz2"
CARD_MEMBER_NAME = "protein_fasta_protein_homolog_model.fasta"
VFDB_URL = "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
BACMET_EXP_URL = "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta"
BACMET_PREDICTED_URL = "http://bacmet.biomedicine.gu.se/download/BacMet2_predicted_database.fasta.gz"
ISFINDER_URL = "https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.fna"
UNIPROT_SPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

_CHUNK_SIZE = 1024 * 1024  # 1MB
# Some sources (confirmed: VFDB's server) return HTTP 403 for urllib's
# default "Python-urllib/x.y" User-Agent, evidently as an anti-scraper
# rule, but allow a normal browser-looking one -- curl works against the
# same URL with its own default UA, and this string reproduces that.
_USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) plasmid-mapper-gen-db-fetch"


class FetchError(Exception):
    """User-facing error with a plain-English message, matching this
    project's convention (pipeline/plasmid_mapper_gen/utils/errors.py's
    PipelineError) of not letting raw tracebacks reach a microbiologist
    running this from the command line.
    """


def _download_to(url: str, dest_path: str, context: str) -> None:
    """Downloads url to dest_path, streaming in chunks (never loading the
    whole response into memory -- UniProt's file alone is ~300MB
    decompressed). Downloads to a temporary sibling file first and only
    renames it to dest_path on success, so a network failure mid-download
    never leaves a truncated file sitting at dest_path looking
    "already downloaded" to a later run's skip-if-present check.
    """
    tmp_fd, tmp_path = tempfile.mkstemp(
        dir=os.path.dirname(dest_path) or ".", prefix=os.path.basename(dest_path) + ".part-"
    )
    try:
        with os.fdopen(tmp_fd, "wb") as tmp_file:
            try:
                request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
                with urllib.request.urlopen(request, timeout=60) as response:
                    shutil.copyfileobj(response, tmp_file, _CHUNK_SIZE)
            except urllib.error.URLError as exc:
                hint = ""
                if url.startswith("http://"):
                    hint = (
                        " This URL is plain HTTP (not HTTPS), which some "
                        "networks/proxies block by policy -- that may be why "
                        "this failed."
                    )
                raise FetchError(f"{context}: failed to download {url}: {exc}.{hint}") from exc
        if os.path.getsize(tmp_path) == 0:
            raise FetchError(f"{context}: downloaded file from {url} is empty")
        os.replace(tmp_path, dest_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _gunzip(gz_path: str, dest_path: str) -> None:
    tmp_fd, tmp_path = tempfile.mkstemp(
        dir=os.path.dirname(dest_path) or ".", prefix=os.path.basename(dest_path) + ".part-"
    )
    try:
        with os.fdopen(tmp_fd, "wb") as tmp_file, gzip.open(gz_path, "rb") as gz_file:
            shutil.copyfileobj(gz_file, tmp_file, _CHUNK_SIZE)
        if os.path.getsize(tmp_path) == 0:
            raise FetchError(f"Decompressing {gz_path} produced an empty file")
        os.replace(tmp_path, dest_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _already_present(dest_path: str, force: bool) -> bool:
    return (not force) and os.path.exists(dest_path) and os.path.getsize(dest_path) > 0


def fetch_card(db_dir: str, force: bool, version: str) -> None:
    dest_path = os.path.join(db_dir, "card.fasta")
    if _already_present(dest_path, force):
        print(f"[card] {dest_path} already present, skipping (use --force to refetch)")
        return

    url = CARD_URL_TEMPLATE.format(version=version)
    print(f"[card] downloading CARD v{version} from {url} ...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        tarball_path = os.path.join(tmp_dir, "card.tar.bz2")
        _download_to(url, tarball_path, "CARD")

        extract_dir = os.path.join(tmp_dir, "extracted")
        os.makedirs(extract_dir, exist_ok=True)
        with tarfile.open(tarball_path, "r:bz2") as tar:
            # filter="data" (safe-extraction filter) was only added as an
            # accepted keyword in Python 3.12 (backported to some patch
            # releases of 3.10/3.11) -- this project's declared minimum is
            # 3.10, so fall back to plain extractall() on older
            # interpreters that reject the kwarg entirely, rather than
            # hard-requiring 3.12+.
            try:
                tar.extractall(extract_dir, filter="data")
            except TypeError:
                tar.extractall(extract_dir)

        member_path = os.path.join(extract_dir, CARD_MEMBER_NAME)
        if not os.path.exists(member_path) or os.path.getsize(member_path) == 0:
            raise FetchError(
                f"CARD v{version}'s archive did not contain the expected "
                f"'{CARD_MEMBER_NAME}' (or it was empty) -- CARD's internal "
                "tarball layout can change between releases; inspect the "
                f"archive at {url} manually to find the right protein "
                "homolog-model FASTA file for this version."
            )
        shutil.copy2(member_path, dest_path)
    print(f"[card] wrote {dest_path}")


def fetch_vfdb(db_dir: str, force: bool) -> None:
    dest_path = os.path.join(db_dir, "vfdb.fasta")
    if _already_present(dest_path, force):
        print(f"[vfdb] {dest_path} already present, skipping (use --force to refetch)")
        return

    print(f"[vfdb] downloading VFDB setA (protein) from {VFDB_URL} ...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        gz_path = os.path.join(tmp_dir, "vfdb.fas.gz")
        _download_to(VFDB_URL, gz_path, "VFDB")
        _gunzip(gz_path, dest_path)
    print(f"[vfdb] wrote {dest_path}")


def fetch_bacmet(db_dir: str, force: bool) -> None:
    dest_path = os.path.join(db_dir, "bacmet.fasta")
    if _already_present(dest_path, force):
        print(f"[bacmet] {dest_path} already present, skipping (use --force to refetch)")
        return

    print(f"[bacmet] downloading experimentally-confirmed set from {BACMET_EXP_URL} ...")
    print(f"[bacmet] downloading predicted set from {BACMET_PREDICTED_URL} ...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        exp_path = os.path.join(tmp_dir, "exp.fasta")
        predicted_gz_path = os.path.join(tmp_dir, "predicted.fasta.gz")
        predicted_path = os.path.join(tmp_dir, "predicted.fasta")
        _download_to(BACMET_EXP_URL, exp_path, "BacMet (experimentally-confirmed)")
        _download_to(BACMET_PREDICTED_URL, predicted_gz_path, "BacMet (predicted)")
        _gunzip(predicted_gz_path, predicted_path)

        # Concatenate both into one bacmet.fasta -- broader recall than
        # either alone; every hit still has to clear the pipeline's own
        # --min-identity/--min-coverage thresholds regardless of which
        # BacMet file it came from, so this doesn't weaken classification
        # quality the way it might if there were no downstream threshold.
        tmp_fd, tmp_combined = tempfile.mkstemp(dir=db_dir, prefix="bacmet.fasta.part-")
        with os.fdopen(tmp_fd, "wb") as combined:
            with open(exp_path, "rb") as f:
                shutil.copyfileobj(f, combined)
            with open(predicted_path, "rb") as f:
                shutil.copyfileobj(f, combined)
        os.replace(tmp_combined, dest_path)
    print(f"[bacmet] wrote {dest_path} (experimentally-confirmed + predicted)")


def fetch_isfinder(db_dir: str, force: bool) -> None:
    dest_path = os.path.join(db_dir, "isfinder.fasta")
    if _already_present(dest_path, force):
        print(f"[isfinder] {dest_path} already present, skipping (use --force to refetch)")
        return

    print(f"[isfinder] downloading nucleotide IS sequences from {ISFINDER_URL} ...")
    _download_to(ISFINDER_URL, dest_path, "ISfinder")
    print(f"[isfinder] wrote {dest_path}")


def fetch_uniprot(db_dir: str, force: bool) -> None:
    dest_path = os.path.join(db_dir, "uniprot_sprot.fasta")
    if _already_present(dest_path, force):
        print(f"[uniprot] {dest_path} already present, skipping (use --force to refetch)")
        return

    print(f"[uniprot] downloading reviewed SwissProt (~300MB decompressed) from {UNIPROT_SPROT_URL} ...")
    print("[uniprot] this is the largest of the 5 databases and may take a while")
    with tempfile.TemporaryDirectory() as tmp_dir:
        gz_path = os.path.join(tmp_dir, "uniprot_sprot.fasta.gz")
        _download_to(UNIPROT_SPROT_URL, gz_path, "UniProt/SwissProt")
        _gunzip(gz_path, dest_path)
    print(f"[uniprot] wrote {dest_path}")


FETCHERS = {
    "card": lambda db_dir, force, args: fetch_card(db_dir, force, args.card_version),
    "isfinder": lambda db_dir, force, args: fetch_isfinder(db_dir, force),
    "vfdb": lambda db_dir, force, args: fetch_vfdb(db_dir, force),
    "bacmet": lambda db_dir, force, args: fetch_bacmet(db_dir, force),
    "uniprot": lambda db_dir, force, args: fetch_uniprot(db_dir, force),
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Downloads the 5 reference databases plasmid-mapper-gen "
            "classifies ORFs against into --db-dir. Safe to re-run: "
            "already-downloaded files are skipped unless --force is given."
        )
    )
    parser.add_argument("--db-dir", required=True, help="Directory to write the database files into")
    parser.add_argument(
        "--force", action="store_true", help="Re-download even if a file already exists at the destination"
    )
    parser.add_argument(
        "--only",
        help=(
            "Comma-separated subset of databases to fetch (choices: "
            + ", ".join(FETCHERS) + "). Default: all 5. "
            "Use this to refresh a single database, e.g. --only uniprot"
        ),
    )
    parser.add_argument(
        "--card-version",
        default="4.0.0",
        help="CARD release version to fetch (default: 4.0.0, pinned rather than 'latest')",
    )
    args = parser.parse_args(argv)

    os.makedirs(args.db_dir, exist_ok=True)

    if args.only:
        requested = [name.strip() for name in args.only.split(",") if name.strip()]
        unknown = [name for name in requested if name not in FETCHERS]
        if unknown:
            print(f"Error: unknown database name(s): {', '.join(unknown)}. Choices: {', '.join(FETCHERS)}", file=sys.stderr)
            return 1
    else:
        requested = list(FETCHERS)

    had_error = False
    for name in requested:
        try:
            FETCHERS[name](args.db_dir, args.force, args)
        except FetchError as exc:
            print(f"Error: {exc}", file=sys.stderr)
            had_error = True

    return 1 if had_error else 0


if __name__ == "__main__":
    sys.exit(main())
