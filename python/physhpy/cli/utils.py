import collections
import csv
import re

Sequence = collections.namedtuple("sequence", ["taxon", "sequence"])


def read_fasta_sequences(filename: str) -> list[Sequence]:
    sequences = {}
    with open(filename, "r") as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(">"):
                taxon = line[1:]
                sequences[taxon] = ""
            else:
                sequences[taxon] += line
    return [Sequence(taxon, sequence) for taxon, sequence in sequences.items()]


def convert_date_to_real(day, month, year):
    if year % 4 == 0:
        days = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    else:
        days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

    for i in range(month - 1):
        day += days[i]

    return (day - 1) / sum(days) + year


def read_dates_from_csv(input_file, date_format=None):
    dates = {}
    with open(input_file) as fp:
        reader = csv.reader(
            fp,
            quotechar='"',
            delimiter=",",
            quoting=csv.QUOTE_ALL,
            skipinitialspace=True,
        )
        for line in reader:
            index_name = line.index("strain")
            index_date = line.index("date")
            break
        for line in reader:
            dates[line[index_name]] = line[index_date]

    if date_format is not None:
        res = re.split(r"[/-]", date_format)
        yy = res.index("yyyy")
        MM = res.index("MM")
        dd = res.index("dd")

        for key, date in dates.items():
            res1 = re.split(r"[/-]", date)
            dates[key] = convert_date_to_real(
                int(res1[dd]), int(res1[MM]), int(res1[yy])
            )
    return dates
