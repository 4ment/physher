"""The readers, against the fluA data in examples/."""

import pytest

from phyforge.io import dates_from_names, decimal_year, read_fasta, taxa_from_tree
from phyforge.io.dates import DateParseError


def test_fasta_and_tree_agree_on_taxa(alignment, tree_file):
    from_alignment = [sequence.taxon for sequence in read_fasta(alignment)]
    from_tree = taxa_from_tree(tree_file)
    assert len(from_alignment) == 69
    assert set(from_alignment) == set(from_tree)


def test_sequences_are_read_whole(alignment):
    sequences = read_fasta(alignment)
    assert all(len(sequence.sequence) == 987 for sequence in sequences)


def test_nexus_extension_is_not_trusted(tree_file):
    # fluA-rooted.nxs is plain newick despite its extension.
    assert taxa_from_tree(tree_file)[0] == "A_Belgium_2_1981"


def test_dates_from_taxon_names(tree_file):
    dates = dates_from_names(taxa_from_tree(tree_file))
    assert min(dates.values()) == 1981.0
    assert max(dates.values()) == 1998.0


def test_missing_date_is_reported():
    with pytest.raises(DateParseError, match="undated"):
        dates_from_names(["undated"])


@pytest.mark.parametrize(
    "day,month,year,expected",
    [(1, 1, 2000, 2000.0), (31, 12, 2001, 2001 + 364 / 365)],
)
def test_decimal_year(day, month, year, expected):
    assert decimal_year(day, month, year) == pytest.approx(expected)


def test_leap_year_has_366_days():
    assert decimal_year(31, 12, 2000) == pytest.approx(2000 + 365 / 366)
