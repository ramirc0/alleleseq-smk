"""Compare a ported output table with the upstream AlleleSeq2 table.

Upstream wrote Python floats with 12 significant digits and R numbers with 15,
so cells are compared as numbers: integers and text exactly, floats within a
relative tolerance.
"""

import math


def read_rows(path, skip_prefix=None):
    """Read a tab-separated table.

    Parameters
    ----------
    path : str
        Table path.
    skip_prefix : str, optional
        Lines starting with this prefix are skipped.

    Returns
    -------
    list of list of str
        Rows, header included.
    """
    with open(path) as handle:
        return [
            line.rstrip("\n").split("\t")
            for line in handle
            if skip_prefix is None or not line.startswith(skip_prefix)
        ]


def cells_match(ours, golden, rtol):
    """Tell whether two cells hold the same value.

    Parameters
    ----------
    ours : str
        Cell of the ported table.
    golden : str
        Cell of the upstream table.
    rtol : float
        Relative tolerance for floats.

    Returns
    -------
    bool
        True when the text matches or both parse to close floats.
    """
    if ours == golden:
        return True
    if "__" in ours and "__" in golden:
        parts_a, parts_b = ours.split("__"), golden.split("__")
        return len(parts_a) == len(parts_b) and all(
            cells_match(a, b, rtol) for a, b in zip(parts_a, parts_b)
        )
    try:
        a, b = float(ours), float(golden)
    except ValueError:
        return False
    if math.isnan(a) and math.isnan(b):
        return True
    return math.isclose(a, b, rel_tol=rtol, abs_tol=1e-300)


def compare_tables(ours_path, golden_path, rtol=1e-11, columns=None, skip_prefix=None):
    """Compare two tables cell by cell.

    Parameters
    ----------
    ours_path : str
        Ported table.
    golden_path : str
        Upstream table.
    rtol : float
        Relative tolerance for floats.
    columns : list of int, optional
        Column indices to compare. Defaults to all.
    skip_prefix : str, optional
        Lines starting with this prefix are skipped in both files.

    Returns
    -------
    list of str
        Human-readable mismatches. Empty when the tables agree.
    """
    ours = read_rows(ours_path, skip_prefix)
    golden = read_rows(golden_path, skip_prefix)
    problems = []
    if len(ours) != len(golden):
        problems.append(f"row count {len(ours)} != {len(golden)}")
    for i, (row_a, row_b) in enumerate(zip(ours, golden)):
        picked = columns if columns is not None else range(max(len(row_a), len(row_b)))
        for j in picked:
            a = row_a[j] if j < len(row_a) else "<missing>"
            b = row_b[j] if j < len(row_b) else "<missing>"
            if not cells_match(a, b, rtol):
                problems.append(f"row {i} col {j}: {a!r} != {b!r}")
                if len(problems) > 20:
                    return problems
    return problems
