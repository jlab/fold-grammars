from math import log, exp
import os
import click
import click_option_group
from hashlib import md5

from parse import *
from execute import *
from tempfile import gettempdir
import pickle

PROGNAME = 'RNAhybrid'
# values taken from RNAhybrid-2.1.2/src/globals.h
DISTRIBUTION = {
    '3utr_fly': {
        'xi_slope': -0.03144,
        'xi_intercept': 0.7201,
        'theta_slope': -0.003429,
        'theta_intercept': 0.01634},
    '3utr_human':{
        'xi_slope': -0.01237,
        'xi_intercept': 1.901,
        'theta_slope': -0.001678,
        'theta_intercept': 0.1349},
    '3utr_worm':{
        'xi_slope': -0.01074,
        'xi_intercept': 1.769,
        'theta_slope': -0.001154,
        'theta_intercept': 0.1419}}

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option(
    '-t', '--target', type=str, required=True, help='An RNA sequence to be searched for miRNA targets. You can either enter a single RNA sequence OR a path to a (multiple) FASTA file.')
@click.option(
    '-q', '--mirna', type=str, required=True, help='A microRNA sequence that shall bind to the target sequence. You can either enter a single RNA sequence OR a path to a (multiple) FASTA file.')
@click_option_group.optgroup.group(
    'p-value calculation', cls=click_option_group.RequiredMutuallyExclusiveOptionGroup)
@click_option_group.optgroup.option(
    '-s', '--set', type=click.Choice(['3utr_fly','3utr_worm','3utr_human']), help='A pre-trained set of organism specific values.')
@click_option_group.optgroup.option(
    '-d', '--distribution', nargs=2, type=float, help='<xi> and <theta> are the position and shape parameters, respectively, of the extreme value distribution assumed for p-value calculation. For example enter "0.91 0.25" without the quotes.')
@click.option(
    '--binPath', type=str, default="", help='%s expects that according Bellman\'s GAP compiled binaries are located in the same directory as the Python wrapper is. Should you moved them into another directory, you must set --binPath to this new location!' % PROGNAME)
@click.option(
    '--filter_minmax_seedlength', type=int, nargs=2, default=(0,9999), help='Minimal and maximal number of base-pairs in the seed helix.')
@click.option(
    '--filter_minmax_mirnabulgelength', type=int, nargs=2, default=(0,9999), help='Minimal and maximal number of bulged bases in microRNA directly after seed helix.')
@click.option(
    '--filter_max_energy', type=int, nargs=1, default=0, help='Maximal free energy (not that stable energies are negative)')
@click.option(
    '--verbose/--no-verbose', type=bool, default=False, help="Print verbose messages.")
@click.option(
    '--cache/--no-cache', type=bool, default=False, help="Cache the last execution and reuse if input does not change. Will store files to '%s'" % gettempdir())
def RNAhybrid(target, mirna, set, distribution, binpath, filter_minmax_seedlength, filter_minmax_mirnabulgelength, filter_max_energy, verbose, cache):
    settings = dict()

    settings['binpath'] = binpath
    settings['program_name'] = PROGNAME

    if filter_minmax_seedlength[0] > filter_minmax_seedlength[1]:
        raise ValueError("minimum for seed length must be smaller than maximum!")
    if filter_minmax_mirnabulgelength[0] > filter_minmax_mirnabulgelength[1]:
        raise ValueError("minimum for bulge must be smaller than maximum!")

    entries_mirnas = []
    if os.path.exists(mirna) and (mirna.upper().replace('A', '').replace('C', '').replace('G', '').replace('U', '') != ""):
        entries_mirnas = list(read_fasta(mirna))
    else:
        entries_mirnas = [('from commandline', mirna)]

    entries_targets = []
    if os.path.exists(target) and (target.upper().replace('A', '').replace('C', '').replace('G', '').replace('U', '') != ""):
        entries_targets = read_fasta(target)
    else:
        entries_targets = [('from commandline', target)]

    # compute maximum duplex energies once per miRNA
    mdes = dict()
    if not distribution and set:
        for entry_mirna in entries_mirnas:
            cmd_mde = compose_call('mde', '', entry_mirna[1], complement(entry_mirna[1]), **settings)
            if verbose:
                print("Binary call: %s" % cmd_mde, file=sys.stderr)
            res_mde = Product(TypeFloat('mde')).parse_lines(execute(cmd_mde))[0]
            mdes[entry_mirna[0]] = res_mde['mde'] / 100

    result_nr = 1
    for entry_target in entries_targets:
        for entry_mirna in entries_mirnas:
            settings['xi'] = 0
            settings['theta'] = 0
            if not distribution and set:
                # estimate evd parameters from maximal duplex energy
                settings['xi'] = DISTRIBUTION[set]['xi_slope'] * mdes[entry_mirna[0]] + DISTRIBUTION[set]['xi_intercept']
                settings['theta'] = DISTRIBUTION[set]['theta_slope'] * mdes[entry_mirna[0]] + DISTRIBUTION[set]['theta_intercept']

            cmd_hybrid = compose_call('khorshid', 'rnahybrid', entry_target[1], entry_mirna[1], **settings)
            if verbose:
                print("Binary call: %s" % cmd_hybrid, file=sys.stderr)
            fp_cache = os.path.join(gettempdir(), md5(cmd_hybrid.encode()).hexdigest() + '.rnahybrid')
            raw_hybrid = []
            if os.path.exists(fp_cache) and cache:
                if verbose:
                    print("Read cached result from file '%s'" % fp_cache, file=sys.stderr)
                with open(fp_cache, 'r') as f:
                    raw_hybrid = f.read().splitlines()
            else:
                raw_hybrid = execute(cmd_hybrid)
                if cache:
                    with open(fp_cache, 'w') as f:
                        f.write('\n'.join(raw_hybrid))
                    if verbose:
                        print("Wrote results into cache file '%s'" % fp_cache, file=sys.stderr)

            #res_stacklen = Product(Product(TypeInt('stacklen'), TypeMFE()), TypeHybrid()).parse_lines(raw_hybrid)
            res_stacklen = Product(Product(TypeKhorshid(), TypeMFE()), TypeHybrid()).parse_lines(raw_hybrid)

            out = sys.stdout
            for answer in sorted(res_stacklen, key=lambda a: (a['stacklen'], a['mfe'])):
                if (answer['stacklen'] < filter_minmax_seedlength[0]) or (answer['stacklen'] > filter_minmax_seedlength[1]):
                    continue
                if (answer['bulgelen'] < filter_minmax_mirnabulgelength[0]) or (answer['bulgelen'] > filter_minmax_mirnabulgelength[1]):
                    continue
                if (answer['mfe'] / 100 > filter_max_energy):
                    continue
                out.write("result no %i\n" % result_nr)
                result_nr += 1
                out.write("target: %s\n" % entry_target[0])
                out.write("length: %i\n" % len(entry_target[1]))
                out.write("miRNA : %s\n" % entry_mirna[0])
                out.write("length: %i\n" % len(entry_mirna[1]))
                out.write("mfe: %.1f kcal/mol\n" % (answer['mfe'] / 100))
                normalized_energy = (answer['mfe'] / 100) / log(len(entry_target[1]) * len(entry_mirna[1]))
                pvalue = 1 - exp(-1 * exp(-1 * ((-1 * normalized_energy - settings['xi']) / settings['theta'])))
                out.write("p-value: %f\n" % pvalue)
                out.write("5' seed length: %s\n" % answer['stacklen'])
                out.write("miRNA bulge length: %s\n" % answer['bulgelen'])
                out.write("\n")
                out.write("position  %i\n" % (answer['target_position']))
                out.write("target 5' " + answer['target_unpaired'] + " 3'\n")
                out.write("          " + answer['target_stacked'] + "\n")
                out.write("          " + answer['pairs'] + "\n")
                out.write("          " + answer['mirna_stacked'] + "\n")
                out.write("miRNA  3' " + answer['mirna_unpaired'] + " 5'\n")
                out.write("\n")
                out.write("\n")


if __name__ == '__main__':
    RNAhybrid()
