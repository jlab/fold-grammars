import sys
from collections import OrderedDict
from typing import Generator, Tuple


def read_fasta(filename:str) -> Generator[Tuple[str, str], None, None]:
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

def read_CT_file(filename:str) -> Generator[Tuple[str, str, dict[int, int]], None, None]:
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

def sort_pairs(pairs: dict[int, int], only_ij=False) -> dict[int, int]:
    sorted_pairs = dict()

    if only_ij:
        # ensure all pairs satisfy i < j
        for i, j in pairs.items():
            if i < j:
                sorted_pairs[i] = j
            else:
                sorted_pairs[j] = i
    else:
        sorted_pairs = pairs
        # ensure that reverse pairs are also included in sorted return set
        sorted_pairs.update({j: i for i,j in sorted_pairs.items()})

    # ensure pairs are sorted by i
    sorted_pairs = OrderedDict({i: j for i, j in sorted(sorted_pairs.items())})

    return sorted_pairs

def get_left_right_positions(pairs: dict[int, int]) -> Tuple[int, int]:
    positions = sorted(list(pairs.keys()) + list(pairs.values()))
    return (positions[0], positions[-1])

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
    orderes_pairs = sort_pairs(pairs, only_ij=True) #OrderedDict({i:j for i,j in pairs.items() if i < j})

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

    if verbose is not None and len(pseudoknotted) > 0:
        print("Warning: %i of %i base-pairs in your structure are pseudoknotted, i.e. crossing." % (len(pseudoknotted), len(nested) + len(pseudoknotted)), file=verbose)

    # return "normal" dictionaries, not ordered ones
    result = {'nested': dict(nested), 'knotted': dict(pseudoknotted)}
    # restore symmetric pairs, e.g. 4: 8 will be complemented by 8: 4
    for (field, bps) in [('nested', nested), ('knotted', pseudoknotted)]:
        for i, j in bps.items():
            result[field][j] = i

    return result

def nested_pairs_to_dotBracket(pairs: dict[int, int], break_positions: [int]=[]) -> str:
    """Takes a dict of pairs and returns a dot-Bracket string.

    Parameters
    ==========
    pairs : dict[int, int]
        A "set" of base-pairs
    break_positions : [int]
        A list of positions that breaks pairs, i.e. affected base-pairs in pairs
        will be represented as . . instead of ( )

    Returns
    =======
        A dot bracket string representation of the nested secondary structure

    Raises
    =====
    ValueError if pairs contain crossing base-pairs
    """

    sorted_pairs = sort_pairs(pairs, only_ij=False)
    dis = disentangle_knots(sorted_pairs)
    if len(dis['knotted']) != 0:
        raise ValueError("Cannot produce dot-Bracket strings for pseudoknotted pair-sets!")

    # determine left and right borders
    (left, right) = get_left_right_positions(sorted_pairs)
    db = ""
    # break_positions might contain partial pairs, we here ensure opening and closing partner are contained
    break_pairs = {i: j for i, j in sorted_pairs.items() if i in break_positions or j in break_positions}
    for x in range(left, right+1, 1):
        y = sorted_pairs.get(x, None)  # get potential pairing partner
        if (y is None) or (x in break_pairs) or (y in break_pairs):
            db += '.'
        elif x < y:
            db += '('
        else:
            db += ')'

    assert db.count('(') == db.count(')')

    return db

def get_minimal_valid_substructure(pairs_full: dict[int, int], pair_subset: dict[int, int]) -> dict[int, int]:
    """Returns all base-pairs of a template structure (pairs_full) that are enclosed by the provided subset.

    Examples
    ========
     pairs_full:    [...]((........)).((((((.....))))))((((..((((((((..(((..((((((((....)))))..)))))).)))))))).((((........))))..............))))((((([...]
     pair_subset:                     ((((((.....))))))((((....(...........................................).................................))))
     result:                          ((((((.....))))))((((..((((((((..(((..((((((((....)))))..)))))).)))))))).((((........))))..............))))

     pairs_full:    [...]((........)).((((((.....))))))((((..((((((((..(((..((((((((....)))))..)))))).)))))))).((((........))))..............))))((((([...]
     pair_subset:                     ((((((.....))))))........(...........................................)
     result:                          ((((((.....))))))((((..((((((((..(((..((((((((....)))))..)))))).)))))))).((((........))))..............))))

    """

    sorted_pairs_full = sort_pairs(pairs_full)
    (left, right) = get_left_right_positions(pair_subset)

    extended_pairs = pair_subset.copy()
    (ext_left, ext_right) = (0, 0)

    while True:
        for i in range(left, right+1, 1):
            if i in pairs_full.keys():
                extended_pairs[i] = pairs_full[i]
                extended_pairs[pairs_full[i]] = i
        (ext_left, ext_right) = get_left_right_positions(extended_pairs)
        if ext_left < left or ext_right > right:
            (left, right) = (ext_left, ext_right)
        else:
            break

    return extended_pairs
