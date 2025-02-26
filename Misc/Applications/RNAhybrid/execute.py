import subprocess
import sys


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



def mainP():
    settings = {'bindir': './x86_64-linux-gnu/'}
    cmd = compose_call('stacklen', 'rnahybrid', 'AUUAAAGGUUUAUACCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUC', 'ACGACAAUCGGGAUCGGGGCGU', **settings)
    run = execute(cmd)


if __name__ == '__main__':
    mainP()
