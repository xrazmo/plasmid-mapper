import importlib.resources

import yaml

from .utils.errors import PipelineError

_DB_KEYS = ["card", "isfinder", "vfdb", "bacmet", "uniprot"]


def _load_yaml_file(path) -> dict:
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise PipelineError(f"{path}: expected a mapping of db name to thresholds")
    return data


def _bundled_defaults() -> dict:
    ref = importlib.resources.files("plasmid_mapper_gen") / "db_thresholds.yaml"
    with importlib.resources.as_file(ref) as path:
        return _load_yaml_file(path)


def load_thresholds(path: str = None) -> dict:
    """Return {db_key: (min_identity, min_coverage)} for all 5 databases.

    Loads the bundled db_thresholds.yaml as the base, then merges a
    user-supplied override file (if `path` is given) on top of it. An
    override file may set only the databases it wants to change --
    missing keys keep the bundled default, and unknown keys are ignored
    (same missing-is-fine philosophy as databases.build_registry(), which
    already treats an absent database file as a valid partial
    configuration rather than an error).
    """
    merged = _bundled_defaults()
    if path is not None:
        override = _load_yaml_file(path)
        for key, value in override.items():
            if key in _DB_KEYS:
                merged[key] = value

    thresholds = {}
    for key in _DB_KEYS:
        entry = merged.get(key, {})
        thresholds[key] = (
            float(entry.get("min_identity", 70.0)),
            float(entry.get("min_coverage", 70.0)),
        )
    return thresholds
