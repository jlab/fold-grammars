import sys

def print_single_answer(target, mirna, answer, number:int, out=sys.stdout):
    out.write("result no %i\n" % number)
    out.write("target: %s\n" % target[0])
    out.write("length: %i\n" % len(target[1]))
    out.write("miRNA : %s\n" % mirna[0])
    out.write("length: %i\n" % len(mirna[1]))
    out.write("mfe: %.1f kcal/mol\n" % (answer['mfe'] / 100))
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

def print_answers(target, mirna, answers, number:int, out=sys.stdout):
    # sort by:
    #   a) ascending p-Value (lower = better)
    #   b) ascending MFE (lower = better, note: stable energies are negative)
    #   c) descending stacklen (larger = better)
    for answer in sorted(answers, key=lambda a: (a['p-value'], a['mfe'], -1*a['stacklen'])):
        inc = print_single_answer(target, mirna, answer, number, out)
        number += inc
    return number

def filter_answers(answers, filter_minmax_seedlength, filter_minmax_mirnabulgelength, filter_max_energy):
    results = []
    for answer in answers:
        if (answer['stacklen'] < filter_minmax_seedlength[0]) or (answer['stacklen'] > filter_minmax_seedlength[1]):
            continue
        if (answer['bulgelen'] < filter_minmax_mirnabulgelength[0]) or (answer['bulgelen'] > filter_minmax_mirnabulgelength[1]):
            continue
        if (answer['mfe'] / 100 > filter_max_energy):
            continue
        results.append(answer)

    return results
