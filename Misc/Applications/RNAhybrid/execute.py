import subprocess
import sys
from collections import OrderedDict


def complement(sequence: str):
    """Returns the complementary sequence."""
    res = ""
    for nuc in sequence.upper():
        res += {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A'}.get(nuc)
    return res


def compose_call(mode: str, grammar: str, inp_target: str, inp_mirna: str, inp_structure: str=None, **kwargs):
    AVAILABLE_MODI = ['stacklen', 'mde', 'khorshid', 'eval']
    AVAILABLE_GRAMMARS = ['rnahybrid', 'microstate']

    cmd = kwargs['binpath']
    cmd += kwargs['program_name']

    if mode not in AVAILABLE_MODI:
        raise ValueError("Unknown mode '%s'! Available modes are '%s'" % (mode, "', '".join(AVAILABLE_MODI)))
    cmd += '_%s' % mode

    if grammar != "":
        if grammar not in AVAILABLE_GRAMMARS:
            raise ValueError("Unknown grammar '%s'! Available gammars are '%s'" % (grammar, "', '".join(AVAILABLE_GRAMMARS)))
        cmd += '_%s' % grammar

    if ('allow_lonelypairs' in kwargs) and (kwargs['allow_lonelypairs'] is True):
        cmd += ' -u 1 '
    cmd += " '%s'" % inp_target.upper().replace('T', 'U')
    if inp_mirna is not None:
        cmd += " '%s'" % inp_mirna.upper().replace('T', 'U')
    if inp_structure is not None:
        cmd += " '%s'" % inp_structure

    return cmd

def execute(cmd):
    with subprocess.Popen([cmd],
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          executable="bash") as call_x:
        if (call_x.wait() != 0):
            out, err = call_x.communicate()
            raise ValueError((
                "SYSTEM CALL FAILED.\n==== STDERR ====\n%s"
                "\n\n==== STDOUT ====\n%s\n") % (
                    err.decode("utf-8", 'backslashreplace'),
                    out.decode("utf-8", 'backslashreplace')))
        else:
            out, err = call_x.communicate()
            return out.decode("utf-8", 'backslashreplace').split('\n')

def read_fasta(filename):
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

def read_CT_file(filename):
    # see https://rna.urmc.rochester.edu/Text/File_Formats.html for format description
    with open(filename, 'r') as f:
        sequence = ""
        structure = dict()
        expected_length = None,
        entry_name = None
        for lnr, line in enumerate(f.readlines()):
            if lnr == 0:
                parts = line.strip().split()
                # 1. Start of first line: number of bases in the sequence
                expected_length = int(parts[0])
                # 2. End of first line: title of the structure
                entry_name = parts[-1]
                continue
            # 3. Each of the following lines provides information about a given base in the sequence. Each base has its own line, with these elements in order:
            #   1 Base number: index n
            #   2 Base (A, C, G, T, U, X)
            #   3 Index n-1
            #   4 Index n+1
            #   5 Number of the base to which n is paired. No pairing is indicated by 0 (zero).
            #   6 Natural numbering. RNAstructure ignores the actual value given in natural numbering, so it is easiest to repeat n here.
            if lnr > expected_length:
                raise ValueError("CT File seems to contain more bases than the %i expected!" % expected_length)
            _, base, _, _, closing_position, opening_position = line.strip().split()
            sequence += base
            if int(closing_position) == 0:
                structure[int(opening_position)] = 'unpaired'
            else:
                structure[int(opening_position)] = int(closing_position)
                structure[int(closing_position)] = int(opening_position)
    if len(structure) != expected_length:
        raise ValueError("Sequence/Structure length not as expected! Please compare to length given in the first line.")

    return {'name': entry_name,
            'sequence': sequence,
            'structure': structure}

def disentangle_knots(structure:dict, verbose=sys.stderr):
    """Disentengles nested and knotted base-pairs"""
    nested = OrderedDict()
    pseudoknotted = dict()
    for (a, b) in sorted(structure.items()):
        if b == 'unpaired':
            continue
        if a > b:
            continue
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

    if len(pseudoknotted) > len(nested):
        nested, pseudoknotted = pseudoknotted, nested

    if verbose is not None:
        print("Warning: %i of %i base-pairs in your structure are pseudoknotted, i.e. crossing." % (len(pseudoknotted), len(nested) + len(pseudoknotted)))

    return {'nested': nested, 'knotted': pseudoknotted}


def mainP():
    settings = {'bindir': './x86_64-linux-gnu/'}
    cmd = compose_call('stacklen', 'rnahybrid', 'AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUC', 'ACGACAAUCGGGAUCGGGGCGU', **settings)
    run = execute(cmd)


if __name__ == '__main__':
    mainP()
