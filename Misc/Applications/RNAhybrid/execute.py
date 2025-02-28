import subprocess
import sys
import os
from tempfile import gettempdir
from hashlib import md5
from output import warning


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
        cmd += ' -u 1'
    cmd += " '%s'" % inp_target.upper().replace('T', 'U')
    if inp_mirna is not None:
        cmd += " '%s'" % inp_mirna.upper().replace('T', 'U')
    if inp_structure is not None:
        cmd += " '%s'" % inp_structure

    return cmd


def execute_subprocess(cmd, verbose=sys.stderr):
    warning("Binary call: %s" % cmd, verbose)
    process = subprocess.Popen([cmd],
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          executable="bash",
                          text=True)
    out = ""
    for oline in iter(process.stdout.readline, ''):
        out += oline
    process.stdout.close()

    err = ""
    for eline in iter(process.stderr.readline, ''):
        err += eline
    process.stderr.close()

    exitCode = process.wait()
    if (exitCode != 0):
        raise ValueError((
                "SYSTEM CALL FAILED.\n==== STDERR ====\n%s"
                "\n\n==== STDOUT ====\n%s\n") % (err, out))
    else:
        return out.split('\n')


def cache_execute(cmd:str, cache, cache_suffix:str='.rnahybrid', verbose=sys.stderr) -> [str]:
    if cache is False:
        return execute_subprocess(cmd, verbose)

    fp_cache = os.path.join(gettempdir(), md5(cmd.encode()).hexdigest() + cache_suffix)
    raw = []
    if os.path.exists(fp_cache):
        warning("Read cached result from file '%s'" % fp_cache, verbose)
        with open(fp_cache, 'r') as f:
            raw = f.read().splitlines()
    else:
        raw = execute_subprocess(cmd, verbose)
        with open(fp_cache, 'w') as f:
            f.write('\n'.join(raw))
        warning("Wrote results into cache file '%s'" % fp_cache, verbose)

    return raw
