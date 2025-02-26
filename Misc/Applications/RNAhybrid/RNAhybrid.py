from math import log, exp
import os
import click
import click_option_group
from hashlib import md5
from multiprocessing import Pool, cpu_count

from parse import *
from execute import *
from tempfile import gettempdir
from output import *
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

def wrap_process(all_data):
    return process_onetarget_onemirna(*all_data)

def process_eval(sequence, dotBracket, verbose, cache, settings):
    # settings for eval binary call
    eval_settings = settings.copy()
    eval_settings.update({'allow_lonelypairs': True})

    cmd_eval = compose_call('eval', 'microstate', sequence, None, inp_structure=dotBracket, **eval_settings)
    fp_cache = os.path.join(gettempdir(), md5(cmd_eval.encode()).hexdigest() + '.rnaeval')
    raw_eval = []
    if os.path.exists(fp_cache) and cache:
        if verbose:
            print("Read cached result from file '%s'" % fp_cache, file=sys.stderr)
        with open(fp_cache, 'r') as f:
            raw_eval = f.read().splitlines()
    else:
        raw_eval = execute(cmd_eval)
        if cache:
            with open(fp_cache, 'w') as f:
                f.write('\n'.join(raw_eval))
            if verbose:
                print("Wrote results into cache file '%s'" % fp_cache, file=sys.stderr)
    res_eval = Product(TypeMFE()).parse_lines(raw_eval)

    return res_eval[0]['mfe']

def process_onetarget_onemirna(entry_target, pos_target, entry_mirna, pos_mirna, mdes, distribution, set, verbose, cache, settings):
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
    # add original positions for each answer of
    #   a) target sequence position in input file
    #   b) mirna sequence position in input file
    #   c) position in gapc raw result
    target_structure = None
    if (len(entry_target) >= 3) and (entry_target[2] is not None):
        # If user provides ONE given secondary structure for target sequence,
        # we get this set of pairs as the third component of the entry_target tuple.
        target_structure = entry_target[2]
        # Pairs are stored as opening: closing base pair positions, i.e. opening < closing.
        # To ease access, we enrich this dict by also adding in closing: opening information.
        target_structure.update({c: o for o, c in target_structure.items()})

    for pos_answer, answer in enumerate(res_stacklen):
        answer['pos_target'] = pos_target
        answer['pos_mirna'] = pos_mirna
        answer['pos_answer'] = pos_answer

        answer['target_name'] = entry_target[0]
        answer['target_length'] = len(entry_target[1])
        answer['mirna_name'] = entry_mirna[0]
        answer['mirna_length'] = len(entry_mirna[1])

        # calculate p-value for answer
        normalized_energy = (answer['mfe'] / 100) / log(answer['target_length'] * answer['mirna_length'])
        answer['p-value'] = 1 - exp(-1 * exp(-1 * ((-1 * normalized_energy - settings['xi']) / settings['theta'])))

        # calculate the free energy that is lost if existing base-pairs in the target structure have to be broken to enable new pairs with miRNA
        if target_structure is not None:
            # get list of target bases involved in novel miRNA binding
            novel_target_pairing_partners = [int(x) for x in answer['listOfTargetPartners'].split(',') if x != ""]
            # obtain original pairing partners in target structure
            original_target_pairs = {o: c for o, c in {x: target_structure.get(x, -1) for x in novel_target_pairing_partners}.items() if c != -1}
            if len(original_target_pairs) > 0:
                # extend set of target pairs to valid sub-structure
                sub_structure = extend_pairs_to_valid_substructure(target_structure, original_target_pairs)
                # sort involved base-pairs for extended sub-structure to enable easy access to extended left and right limits
                h = sorted(list(sub_structure.keys()) + list(sub_structure.values()))
                # get a valid short dot-Bracket representation of that part of the user given target secondary structure that is involved in miRNA binding
                sub_structure_dotBracket = get_sub_dotBracket_structure(target_structure, h[0], h[-1])
                assert sub_structure_dotBracket.count('(') == sub_structure_dotBracket.count(')')
                sub_sequence = entry_target[1][h[0]-1:h[-1]]
                answer['original_target_substructure_energy'] = process_eval(sub_sequence, sub_structure_dotBracket, verbose, cache, settings)

                # same as "sub_structure_dotBracket", but those base-pairs that are involved in miRNA binding are broken, i.e. converted to .
                broken_sub_structure_dotBracket = get_sub_dotBracket_structure(target_structure, h[0], h[-1], novel_target_pairing_partners)
                assert broken_sub_structure_dotBracket.count('(') == broken_sub_structure_dotBracket.count(')')
                answer['broken_target_substructure_energy'] = process_eval(sub_sequence, broken_sub_structure_dotBracket, verbose, cache, settings)
                # if (sub_sequence, broken_sub_structure_dotBracket) not in evaluated_structures:
                #     cmd_eval = compose_call('eval', 'microstate', sub_sequence, None, inp_structure=broken_sub_structure_dotBracket, **eval_settings)
                #     raw_eval = execute(cmd_eval)
                #     res_eval = Product(TypeMFE()).parse_lines(raw_eval)
                #     evaluated_structures[(sub_sequence, broken_sub_structure_dotBracket)] = res_eval[0]
                # answer['broken_target_substructure_energy'] = evaluated_structures[(sub_sequence, broken_sub_structure_dotBracket)]['mfe']

    return res_stacklen

def extend_pairs_to_valid_substructure(full_structure, pair_subset):
    """Recursively check if some of the pairs might extend into regions not yet within the sub-structure limits"""

    curr_list_of_paired_bases = sorted(list(pair_subset.keys()) + list(pair_subset.values()))
    orig_left, orig_right = curr_list_of_paired_bases[0], curr_list_of_paired_bases[-1]

    curr_left, curr_right = orig_left, orig_right
    for i in range(orig_left, orig_right+1, 1):
         if i in full_structure.keys():
             j = full_structure[i]
             pair_subset[i] = j
             pair_subset[j] = i
             if j < curr_left:
                 curr_left = j
             if j > curr_right:
                 curr_right = j

    if (orig_left == curr_left) and (orig_right == curr_right):
        return pair_subset
    else:
        for i in range(curr_left, orig_left+1, 1):
            if i in full_structure.keys():
                pair_subset[i] = full_structure[i]
        for i in range(orig_right, curr_right+1, 1):
            if i in full_structure.keys():
                pair_subset[full_structure[i]] = i

        return extend_pairs_to_valid_substructure(full_structure, pair_subset)

def get_sub_dotBracket_structure(structure, left_border:int, right_border:int, break_pairs=[]):
    """Construct a dot-Bracket string for a part of the base-pairs given in the "structure" dictionary,
       but only for the interval between left_border and right_border.
       Do NOT add a () pair, if one of the partners is in the break_pairs list.
    """
    dot = ""
    for i in range(left_border, right_border+1, 1):
        if i not in structure.keys():
            dot += "."
        else:
            j = structure[i]
            if (j < left_border) or (i in break_pairs) or (j in break_pairs):
                dot += '.'
            else:
                if i < j:
                    dot += '('
                else:
                    dot += ')'
    return dot

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click_option_group.optgroup.group(
    'Target Input', cls=click_option_group.RequiredMutuallyExclusiveOptionGroup)
@click_option_group.optgroup.option(
    '-t', '--target', type=str, help='An plain (i.e. no header/name) RNA sequence to be searched for miRNA targets.')
@click_option_group.optgroup.option(
    '-tf', '--target_file', type=click.Path(exists=True, dir_okay=False), help='Read one or multiple target RNA sequences from a FASTA file.')
@click_option_group.optgroup.option(
    '-tct', '--target_CT_file', type=click.Path(exists=True, dir_okay=False), help='Read a target sequence/structure pair from a CT file and relate hybrid energies to energy necessary to break given target structure.')
@click_option_group.optgroup.group(
    'miRNA Input', cls=click_option_group.RequiredMutuallyExclusiveOptionGroup)
@click_option_group.optgroup.option(
    '-q', '--mirna', type=str, help='A plain (i.e. no header/name) microRNA sequence that shall bind to a target sequence.')
@click_option_group.optgroup.option(
    '-qf', '--mirna_file', type=click.Path(exists=True, dir_okay=False), help='Read one or multiple miRNA sequences from a FASTA file.')
@click_option_group.optgroup.group(
    'p-value calculation', cls=click_option_group.RequiredMutuallyExclusiveOptionGroup)
@click_option_group.optgroup.option(
    '-s', '--set', type=click.Choice(['3utr_fly','3utr_worm','3utr_human']), help='A pre-trained set of organism specific values.')
@click_option_group.optgroup.option(
    '-d', '--distribution', nargs=2, type=float, help='<xi> and <theta> are the position and shape parameters, respectively, of the extreme value distribution assumed for p-value calculation. For example enter "0.91 0.25" without the quotes.')
@click.option(
    '--binPath', type=str, default="", help='%s expects that according Bellman\'s GAP compiled binaries are located in the same directory as the Python wrapper is. Should you moved them into another directory, you must set --binPath to this new location!' % PROGNAME)
@click.option(
    '--filter-minmax-seedlength', type=int, nargs=2, default=(0,9999), help='Minimal and maximal number of base-pairs in the seed helix.')
@click.option(
    '--filter-minmax-mirnabulgelength', type=int, nargs=2, default=(0,9999), help='Minimal and maximal number of bulged bases in microRNA directly after seed helix.')
@click.option(
    '--filter-max-energy', type=int, nargs=1, default=0, help='Maximal free energy (not that stable energies are negative)')
@click.option(
    '--filter-max-pvalue', type=click.FloatRange(0, 1), default=1, help='Maximal p-value (the smaller the less likely the result is due to rare chance)')
@click.option(
    '--verbose/--no-verbose', type=bool, default=False, help="Print verbose messages.")
@click.option(
    '--cache/--no-cache', type=bool, default=False, help="Cache the last execution and reuse if input does not change. Will store files to '%s'" % gettempdir())
@click.option(
    '--stream-output/--no-stream-output', type=bool, default=False, help="Internally, computation is done per target per miRNA. By default, results are reported once ALL computations are done. With this setting you can make %s print results as soon as they are available - with the downside that ordering might be a bit more unintuitive." % (PROGNAME))
@click.option(
    "--num-cpus", type=click.IntRange(1, cpu_count()), default=1, help="Number of CPU-cores to use. Default is 1, i.e. single-threaded. Note that --stream-output cannot be used as concurrent sub-tasks would overwrite their results.")
@click.option(
    '--sam', type=click.File('w'), required=False, help="Provide a filename if you want to make %s report results as a *.sam file, such that you can view hit positions in a genome browser like IGV." % PROGNAME)
def RNAhybrid(target, target_file, target_ct_file, mirna, mirna_file, set, distribution, binpath, filter_minmax_seedlength, filter_minmax_mirnabulgelength, filter_max_energy, filter_max_pvalue, verbose, cache, stream_output, num_cpus, sam):
    settings = dict()

    settings['binpath'] = binpath
    settings['program_name'] = PROGNAME

    if filter_minmax_seedlength[0] > filter_minmax_seedlength[1]:
        raise ValueError("minimum for seed length must be smaller than maximum!")
    if filter_minmax_mirnabulgelength[0] > filter_minmax_mirnabulgelength[1]:
        raise ValueError("minimum for bulge must be smaller than maximum!")

    entries_mirnas = []
    if mirna_file:
        entries_mirnas = list(read_fasta(mirna_file))
    else:
        entries_mirnas = [('from commandline', mirna)]

    entries_targets = []
    if target_file:
        entries_targets = read_fasta(target_file)
    elif target_ct_file:
        trgt = read_CT_file(target_ct_file)
        entries_targets = [(
            trgt['name'], trgt['sequence'],
            disentangle_knots(trgt['structure'])['nested'])]
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
    answers = []
    tasks = []
    targets = []
    total_mirnas = 0
    for pos_target, entry_target in enumerate(entries_targets):
        targets.append([entry_target[0], len(entry_target[1])])
        for pos_mirna, entry_mirna in enumerate(entries_mirnas):
            if pos_target == 0:
                total_mirnas += 1
            args = (entry_target, pos_target, entry_mirna, pos_mirna, mdes, distribution, set, verbose, cache, settings)
            if num_cpus == 1:
                res_stacklen = process_onetarget_onemirna(*args)
                if stream_output:
                    answers = res_stacklen
                    answers = filter_answers(answers, filter_minmax_seedlength, filter_minmax_mirnabulgelength, filter_max_energy, filter_max_pvalue)
                    result_nr += print_answers(answers, result_nr)
                else:
                    answers.extend(res_stacklen)
            else:
                tasks.append(args)

    if tasks != []:
        pool = Pool(num_cpus)
        for res in zip(pool.map(wrap_process, tasks)):
            answers.extend(res[0])

    if (not stream_output) or (tasks != []):
        answers = filter_answers(answers, filter_minmax_seedlength, filter_minmax_mirnabulgelength, filter_max_energy, filter_max_pvalue)
        result_nr += print_answers(answers, result_nr, out=sys.stdout)

    print_summary(answers, len(targets), total_mirnas)
    if sam:
        print_sam(answers, targets, sam)

if __name__ == '__main__':
    RNAhybrid()
