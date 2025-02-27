import sys
from pathlib import Path
from io import TextIOWrapper, StringIO

def print_single_answer(answer, number:int, out=sys.stdout):
    out.write("result no %i\n" % number)
    out.write("target: %s\n" % answer['target_name'])
    out.write("length: %i\n" % answer['target_length'])
    out.write("miRNA : %s\n" % answer['mirna_name'])
    out.write("length: %i\n" % answer['mirna_length'])
    out.write("mfe: %.1f kcal/mol\n" % (answer['mfe'] / 100))
    if 'broken_target_substructure_energy' in answer.keys():
        out.write('energy lost to make target accessible: %.1f kcal/mol\n' % ((answer['broken_target_substructure_energy'] - answer['original_target_substructure_energy']) / 100))
    out.write("p-value: %f\n" % answer["p-value"])
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

    return 1

def print_answers(answers, number:int, out=sys.stdout):
    # sort by:
    #   a) ascending p-Value (lower = better)
    #   b) ascending MFE (lower = better, note: stable energies are negative)
    #   c) descending stacklen (larger = better)
    for answer in sorted(answers, key=lambda a: (a['p-value'], a['mfe'], -1*a['stacklen'])):
        inc = print_single_answer(answer, number, out)
        number += inc
    return number

def filter_answers(answers, filter_minmax_seedlength, filter_minmax_mirnabulgelength, filter_max_energy, filter_max_pvalue):
    results = []
    for answer in answers:
        if (answer['stacklen'] < filter_minmax_seedlength[0]) or (answer['stacklen'] > filter_minmax_seedlength[1]):
            continue
        if (answer['bulgelen'] < filter_minmax_mirnabulgelength[0]) or (answer['bulgelen'] > filter_minmax_mirnabulgelength[1]):
            continue
        if (answer['mfe'] / 100 > filter_max_energy):
            continue
        if (answer['p-value'] > filter_max_pvalue):
            continue
        results.append(answer)

    return results

def print_summary(answers, total_targets, total_mirnas, out=sys.stdout):
    out.write('======= Summary =======\n')
    out.write('%i hits reported in total.\n' % len(answers))

    targets = dict()
    order_targets = dict()  # original order in input file
    mirnas = dict()
    order_mirnas = dict()
    for a in answers:
        for (field, d, o) in zip(['target', 'mirna'], [targets, mirnas], [order_targets, order_mirnas]):
            if a['%s_name' % field] not in d:
                d[a['%s_name' % field]] = 0
                o[a['%s_name' % field]] = a['pos_%s' % field]
            d[a['%s_name' % field]] += 1

    for (field, d, o, t) in zip(['target', 'mirna'], [targets, mirnas], [order_targets, order_mirnas], [total_targets, total_mirnas]):
        out.write('Number hits per %s (%i %s with zero hits):\n' % (field, t - len(d), field))
        for name, pos in sorted(o.items(), key=lambda item: item[1]):
            out.write('  %i hits for "%s"\n' % (d[name], name))

def cigar(answer):
    long_cigar = ""
    for i in range(len(answer['target_unpaired'])):
        if (answer['pairs'][i] != " "):
            long_cigar += '='
        else:
            if answer['mirna_stacked'][i] == '-':
                long_cigar += 'D'
            elif answer['mirna_stacked'][i] != ' ':
                long_cigar += 'I'
            else:
                if answer['mirna_unpaired'][i] == '-':
                    long_cigar += 'D'
                elif answer['mirna_unpaired'][i] != ' ':
                    long_cigar += 'X'
                else:
                    long_cigar += 'N'

    lastchar = ""
    cigar = ""
    counts = 0
    for i in range(len(long_cigar)):
        if long_cigar[i] == lastchar:
            counts += 1
        else:
            if lastchar != "":
                cigar += '%i%s' % (counts, lastchar)
            lastchar = long_cigar[i]
            counts = 1
    if lastchar != "":
        cigar += '%i%s' % (counts, lastchar)

    return cigar


def hit_2_sam(answer):
    mirna_sequence = ""
    for i in range(len(answer['mirna_unpaired'])):
        if answer['mirna_unpaired'][i] == " ":
            mirna_sequence += answer['mirna_stacked'][i]
        else:
            mirna_sequence += answer['mirna_unpaired'][i]
    mirna_sequence = mirna_sequence.replace('-', '').strip()

    sam = []
    sam.append('%s:%i' % (answer['mirna_name'], answer['pos_answer']))  # QNAME
    sam.append("0")  # FLAG
    sam.append(answer['target_name'])  # RNAME
    sam.append(str(answer['target_position']))  # POS
    sam.append(str((2**8)-1))  # MAPQ
    #sam.append('%iM' % len(target_sequence)) #hit['cigar'])  # CIGAR
    sam.append(cigar(answer)) #hit['cigar'])  # CIGAR
    sam.append('*')  # RNEXT
    sam.append(str(0))  # PNEXT
    sam.append(str(0))  # TLEN
    sam.append(mirna_sequence)  # SEQ
    sam.append('*')  # QUAL
    sam.append('RG:Z:RNAhybrid')
    sam.append('FE:f:%f' % (answer['mfe'] / 100))
    sam.append('PV:f:%f' % answer['p-value'])
    sam.append('SL:i:%i' % answer['stacklen'])
    sam.append('BL:i:%i' % answer['bulgelen'])
    sam.append("xa:Z:target 5' " + answer['target_unpaired'] + " 3'")
    sam.append("xb:Z:          " + answer['target_stacked'])
    sam.append("xc:Z:          " + answer['pairs'])
    sam.append("xd:Z:          " + answer['mirna_stacked'])
    sam.append("xe:Z:miRNA  3' " + answer['mirna_unpaired'] + " 5'")

    return ('\t'.join(sam))

def print_sam(answers, targets, out=sys.stdout):
    out.write('\t'.join(["@HD", "VN:1.0", "SO:unsorted"])+"\n")

    for (name, length) in targets:
        out.write('\t'.join(["@SQ", "SN:%s" % name, "LN:%i" % length]) + '\n')
    for answer in answers:
        out.write(hit_2_sam(answer) + "\n")

def warning(msg, target, linebreak=True):
    if target is not None:
        if isinstance(target, Path):
            with open(target, 'a') as f:
                f.write(msg)
                if linebreak:
                    f.write('\n')
        elif isinstance(target, TextIOWrapper) or isinstance(target, StringIO):
            print(msg, file=target, end="")
            if linebreak:
                print("", file=target)
