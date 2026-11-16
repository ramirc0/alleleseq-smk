"""Sheet validation in common.smk, checked by a dry run of the integration workflow.

Snakemake lives in the launcher env, not the test env, so the dry run goes
through ``pixi run`` on the root manifest.
"""

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / ".tests/integration/config"


def dry_run(**sheets):
    overrides = [f"{key}={path}" for key, path in sheets.items()]
    return subprocess.run(
        [
            "pixi",
            "run",
            "--manifest-path",
            str(ROOT / "pixi.toml"),
            "snakemake",
            "--snakefile",
            str(ROOT / "workflow/Snakefile"),
            "--directory",
            str(ROOT / ".tests/integration"),
            "--dry-run",
            "--quiet",
            "--config",
            *overrides,
        ],
        capture_output=True,
        check=False,
        text=True,
    )


def duplicated_first_row(sheet, tmp_path):
    header, first, *rest = (CONFIG / sheet).read_text().splitlines()
    path = tmp_path / sheet
    path.write_text("\n".join([header, first, first, *rest]) + "\n")
    return path


def test_integration_sheets_pass():
    result = dry_run()
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    ("key", "sheet", "column"),
    [
        ("individuals", "individuals.tsv", "individual_id"),
        ("samples", "samples.tsv", "sample_id"),
    ],
)
def test_duplicate_ids_are_rejected(key, sheet, column, tmp_path):
    result = dry_run(**{key: duplicated_first_row(sheet, tmp_path)})
    assert result.returncode != 0
    assert f"duplicate {column}" in result.stdout + result.stderr
