import sys
from collections import OrderedDict
from typing import Generator


def read_fasta(filename:str) -> Generator[(str, str)]:
    """Read a (multiple) fastA file"""
    with open(filename, 'r') as f:
        sequence = ""
        header = None
        for line in f.readlines():
            if line.startswith(';'):
                continue
            if line.strip() == "":
                continue
            if line.startswith('>'):
                if (header is not None) and (len(sequence) > 0):
                    res = (header, sequence)
                    sequence = ""
                    yield res
                header = line[1:].strip()
            else:
                sequence += line.strip()
        if (header is not None) and (len(sequence) > 0):
            yield (header, sequence)

def read_CT_file(filename:str) -> Generator[(str, str, dict[int, int])]:
    """Read a (multiple) CT (Connectivity Table) file.

    Parameters
    ----------
    filename : str
        Filepath to input file.

    Returns
    -------
    A generator of triples, where components are
        1: a header descriptor
        2: the nucleotide sequence
        3: all base pairs. Note, each pair i to j is contained as i: j AND j: i
           Further note that number of elements are NOT equal to sequence length,
           as unpaired bases are NOT represented
    """
    # see https://rna.urmc.rochester.edu/Text/File_Formats.html for format description
    with open(filename, 'r') as f:
        next_lines_till_entry_end = 0
        header = None
        sequence = None
        pairs = None
        expected_structure_length = None
        for line in f.readlines():
            if next_lines_till_entry_end == 0:
                # yield current entry
                if header is not None:
                    yield (header, sequence, pairs)

                # new entry
                parts = line.strip().split()
                expected_structure_length = int(parts[0])
                header = ' '.join(parts[1:])
                sequence = ""
                pairs = dict()
                next_lines_till_entry_end = expected_structure_length
            else:
                # 3. Each of the following lines provides information about a given base in the sequence. Each base has its own line, with these elements in order:
                #   1 Base number: index n
                #   2 Base (A, C, G, T, U, X)
                #   3 Index n-1
                #   4 Index n+1
                #   5 Number of the base to which n is paired. No pairing is indicated by 0 (zero).
                #   6 Natural numbering. RNAstructure ignores the actual value given in natural numbering, so it is easiest to repeat n here.
                parts = line.strip().split()

                sequence += parts[1]  # the nucleotide
                curr_base_position = int(parts[0])  # the current nucleotide position
                partner_base_position = int(parts[-2])  # the nucleotide position this nucleotide is paired to. 0 if unpaired
                if partner_base_position > 0:
                    # do not explicitely store unpaired positions
                    pairs[curr_base_position] = partner_base_position

                next_lines_till_entry_end -= 1
        if next_lines_till_entry_end == 0:
            if header is not None:
                yield (header, sequence, pairs)
        else:
            if expected_structure_length - next_lines_till_entry_end > 0:
                raise ValueError("One of the entries in your CT file '%s' seems to miss %i of %i lines." % (filename, next_lines_till_entry_end, expected_structure_length))
            if expected_structure_length == next_lines_till_entry_end:
                raise ValueError("One of the entries in your CT file '%s' seems to have surplus lines. Carefully check the header of the entries and compare to sequence lengths!")

def disentangle_knots(pairs:dict[int, int], verbose=sys.stderr) -> dict[str, dict[int, int]]:
    """Remove minimal number of base-pairs to ensure remaining set is non-crossing.

    Parameters
    ----------
    pairs : dict
        The base-pairs of a secondary structure.
    verbose : stream
        Stream to which warnings are written.

    Returns
    -------
    A dictionary with two keys:
        'nested': dict(int, int) all purely nested base-pairs, including symmetric ones
        'knotted': dict(int, int) those base-pairs, identified as crossing, inclusing symmetric ones
    """
    # only take pairs where i < j is satisfied + order ascendingly
    orderes_pairs = OrderedDict({i:j for i,j in pairs.items() if i < j})

    nested = OrderedDict()
    pseudoknotted = OrderedDict()
    for (a, b) in orderes_pairs.items():
        crossing = False
        for (c, d) in nested.items():
            if d < a:  # 1. we know that c < d AND we know that 5' pairs are all nested
                continue
            if c < a < d < b:
                crossing = True
                pseudoknotted[a] = b
                break
        if crossing is False:
            nested[a] = b

    # flip sets if more pairs occure in "pseuoknotted"
    if len(pseudoknotted) > len(nested):
        nested, pseudoknotted = pseudoknotted, nested

    if verbose is not None:
        print("Warning: %i of %i base-pairs in your structure are pseudoknotted, i.e. crossing." % (len(pseudoknotted), len(nested) + len(pseudoknotted)))

    # return "normal" dictionaries, not ordered ones
    result = {'nested': dict(nested), 'knotted': dict(pseudoknotted)}
    # restore symmetric pairs, e.g. 4: 8 will be complemented by 8: 4
    for (field, bps) in [('nested', nested), ('knotted', pseudoknotted)]:
        for i, j in bps.items():
            result[field][j] = i

    return result
