"""Read SNV and region count tables for the statistics steps.

SNV tables start with ``#chr`` and hold base counts ``cA`` to ``cN``.
Region tables start with ``#region`` and hold ``hap1_count`` and
``hap2_count``.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class CountTable:
    """A counts table split into its header and rows.

    Attributes
    ----------
    header : list of str
        Column names.
    rows : list of list of str
        Data rows in file order.
    """

    header: list
    rows: list

    @property
    def is_region(self):
        """Tell whether the table is region-level.

        Returns
        -------
        bool
            True for region tables, False for SNV tables.
        """
        return self.header[0] == "#region"

    def column(self, name):
        """Return one column as strings.

        Parameters
        ----------
        name : str
            Column name.

        Returns
        -------
        list of str
            The column's cells.
        """
        index = self.header.index(name)
        return [row[index] for row in self.rows]


def read_count_table(path):
    """Read a tab-separated counts table.

    Parameters
    ----------
    path : str
        Table path.

    Returns
    -------
    CountTable
        Header and rows.
    """
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in handle]
    return CountTable(header, rows)


def two_largest(table):
    """Return the largest and second-largest allele count of each row.

    Parameters
    ----------
    table : CountTable
        SNV rows use ``cA`` to ``cT``. Region rows use the two haplotypes.

    Returns
    -------
    higher : numpy.ndarray of int
        Largest count per row.
    lower : numpy.ndarray of int
        Second-largest count per row.
    """
    names = (
        ["hap1_count", "hap2_count"] if table.is_region else ["cA", "cC", "cG", "cT"]
    )
    counts = np.array([table.column(n) for n in names], dtype=np.int64).T.reshape(
        -1, len(names)
    )
    ordered = np.sort(counts, axis=1)
    return ordered[:, -1], ordered[:, -2]


def observed_pvalues(table):
    """Return the binomial p-value column.

    Parameters
    ----------
    table : CountTable
        SNV or region table.

    Returns
    -------
    numpy.ndarray of float
        P-values in row order.
    """
    return np.array(table.column("p_binom"), dtype=np.float64)


def totals_and_ratios(table):
    """Return each row's allele total and tested allele ratio.

    Parameters
    ----------
    table : CountTable
        SNV rows use ``sum_ref_n_alt_cnts`` and ``ref_allele_ratio``. Region
        rows use the haplotype sum and ``hap1_allele_ratio``.

    Returns
    -------
    totals : numpy.ndarray of int
        Reads on the two alleles.
    ratios : numpy.ndarray of float
        Reference (SNV) or hap1 (region) allele ratio.
    """
    if table.is_region:
        totals = np.array(table.column("hap1_count"), dtype=np.int64) + np.array(
            table.column("hap2_count"), dtype=np.int64
        )
        ratios = np.array(table.column("hap1_allele_ratio"), dtype=np.float64)
    else:
        totals = np.array(table.column("sum_ref_n_alt_cnts"), dtype=np.int64)
        ratios = np.array(table.column("ref_allele_ratio"), dtype=np.float64)
    return totals, ratios
