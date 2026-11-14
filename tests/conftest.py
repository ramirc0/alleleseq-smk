"""Shared fixtures for the unit tests."""

import pytest


@pytest.fixture
def write(tmp_path):
    """Return a helper that writes text to a file under tmp_path."""

    def _write(name, text):
        path = tmp_path / name
        path.write_text(text)
        return str(path)

    return _write
