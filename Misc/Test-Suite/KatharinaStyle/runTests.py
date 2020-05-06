#!/homes/kmaibach/miniconda3/bin/python

import unittest
import os
import subprocess
import tempfile
from collections import OrderedDict
import shutil
import sys


class Test(unittest.TestCase):

    def setUp(self):
        self.fp_workdir = tempfile.mkdtemp()

        # Ordered dictionary with 2nd structure counts ans mfe-values for every sequence
        self.sequences = OrderedDict([
            ("CCCAA+AAAUGGG", ["9", "-410"]),
            ("CCCCCCA+AUCCCCCAAAGGGGGGGGGGGG", ["22979", "-2140"]),
            ("CCC+GGG", ["6", "-250"]),
            ("UCC+GGG", ["6", "-20"]),
            ("CCU+GGG", ["6", "-80"]),
            ("CCC+ACCCAAAGGGGGG", ["66", "-370"]),
            ("CCU+ACCCAAAGGGGGG", ["68", "-260"]),
            ("CCC+AUCCAAAGGGGGG", ["66", "-250"]),
            ("CCC+ACCCCCCAAAGGGGGGGGG", ["1196", "-1360"]),
            ("CCC+AUCCCCCAAAGGGGGGGGG", ["1196", "-1150"]),
            ("CCCCCCCCCAAAGGGGGGA+GGG", ["1196", "-1360"]),
            ("CCCUCCCCCAAAGGGGGGA+GGG", ["1378", "-1340"]),
            ("CCCAACCCCCCAAAGGGGGGA+GGG", ["816", "-1360"]),
            ("CCCAAUCCCCCAAAGGGGGGA+GGG", ["990", "-1340"]),
            ("CCCCAACCCCAA+AGGGGAAGGGGA", ["307", "-1440"]),
            ("CCACCAAAGGACCAAAGGA+GG", ["19", "0"]),
            ("GCAUCUGCGUAUCUCGGAUAAUGCAAC+UUGAUAACUAUGC", ["29096", "-540"]),
            ("CAAUUUAAUCCCCACCGGUCACU+AAAGCCAUUAGG", ["1404", "-220"]),
            # use only for test 2
            ("GCAUCUGCGUAUCUCGGAUAAUGCAAC+UUGAUAACUAUGCUAUGG", ["177456", "-540"]),
            ("GGUUUGGGUCUUAAAGGGGCUCUGGGAUACUCCCA+CGGACC", ["479122", "-1170"]),
            ("GGUUUAUCGGGUCUUAAAGGGGCUCUGGGAUACUCCCA+CGGACC", ["1026645", "-1170"]),
            ("GGUUUAUCGGGUCUUAAAGGGGCUCUGGGAUACAUAUUCCCA+CGGACC", ["2242901", "-1280"]),
            ("GGUUUAUCGGGUGCUUAAAGGGGGCUCUGGGAUACAUAUUCCCA+CGGACC", ["4219244", "-1270"]),
            ("GGUUUAUUAGGAUGUCGGGUGCUUAAAGGGGGCUCACGUCUUGGGAUACAUAUUCCCA+CGGACC", ["660278945", "-2000"]),
            ("GGACAGACUACAUAGCGAUAUUCGCGAGGGACGGUCCAAGCCAUAAAGAAGACAUUAUGCCAACG+AUCGGUAGAGGCGCUUGCUGAGAUAUCUUUGGUCUUCUAGCGCCUAAGGGCUUCCACACUCCUUGCAUCGUAAAAAUCUAAGAA", ["12412333775906063261", "-3300"]),
            ("UAGGAAGAUCGAGUGGCGGCUUCAACCAGGUCGACCCUGUAAAGGGAAGCUUAAAAUCAAGGACCUACCG+UAAUCGUAACUACGCACUAAGCGAAUACGUGGCGAGCCGGA", ["8615975181702", "-2130"]),
            ("GGCCACCCUUGACCGUUUUCCCGCUAGCCAACUUACGUGGUGGGAGGACGGUCAACUAUGUGGACAACUAUAUAGUGGAGCGUAGCGACGGUCUGUUGGUACCAGUGU+CCAACAGGCUGUCCGCUCCUGGCC", ["7930496389100622030", "-7080"]),
            ("GGCCAUGACCGUUUUCCCGCUAGCCAACUUACGUGGUGGGAGGACGGUCAACUAUGUGGACAACUAUAUAGUGGAGCGGACGGUCUGUUGGUACCAGUGU+CCAACAGGCUGUCCGCUCCUGGCC", ["560970781638806571", "-7770"]),
            ("UGCUUUGGUCUCUUAGGUUACAGACCAGAUAACGUGCUCUCUUCUCAAGGUGUCAUGGAAAAACGCUUACUAAGAUGAGGAGCGUUGCAUCUGCGUAUCUCGGAUAAUGCAAC+UUGAUAACUAUGCUAUGG", ["3087043724457592213", "-2510"]),
            ("GGUUUAUUAGGAUGUCGGGUGCUUAAAGGGGGCUCACGUCUAGCACGCUGGGAUACAUAUUCCCAAUUGUU+CGGACC", ["49821805198", "-2120"]),
            ("GGUUUAUUACGAAAGGAUGUCGGGUGCUUAAAGGGGGCUCACGUCUGACGAGCACGCUGGGAUACAUAUUCCCAAUUGUU+CGGACC", ["766440248432", "-2190"]),
            ("GGUUUAUUACGAAAGGAUGUCGGGUGCUUAAAGGGGGCUCACGUCUGACGCCACGCUCAUGGAGCACGCUGGGAUACAUAUUCCCAAUUGUUAGGGCCUGGAAGAC+CGGACC", ["3211827845108803", "-2910"]),
            ("UGGUCUCGGUUUAUUACGAAAGGAUGUCGGGUGCUUAAAGGGGGCUCACGUCUGACGCCACGCUCAUGGAGCACGCUGGGAUACAUAUUCCCAAUUGUUAGGGCCUGGAAGAC+CGGACCGACCAAGCAU", ["709508811230949176", "-3580"]),
            ("GUAGUUGACAGCCCUGAAAUUCGGUCG+CUUUUAAUUAACUCCUUCCUCUUCGGACCGAUGUGCCCUUGGAUAUCCUAUAUGGAAUUUGGUGUCAGGAUUGUGCCAAGGCCACGAACACAGGAAUGCCAUACAA", ["462656285088435985", "-3150"]),
            ## ("CCACGGAGUUGUACAUCUCCUGUGCGAUGCCCCGCAUCCUGAAUAGGGCAGGUUGCCAGUGGAGGGUUGGGGAUACUUCCGAA+UCCAUAUCCACUCGUGACCGGAGGUUAGUACUUUACGAGCGAAAGGGUGCCCGCGCUUCUUUUCCGGACCCCGCACCACGCAUGUAGAGUACUUGGGAACGAAAGGCCGUAACACGGUGGGACAAAAGCAGCUGCCUAUACGCCGGUUUACCAAGGUCCCCUGGUUGAAUAUCAUGUCUGUGUUGAA", 223925561507966359),
            ## ("GGCACACCCUUGACCGUUUCCUCGUUCCGCUAGCCAACUUACGUGGUGGCAGAUUGGACGGUCAACUAUAUGUGGACAACUAUAAAUUAGUGGAGCGUAGCGACGGAUCUGUAUUGGUACCAGUGU+CCAAACCAGGCUGUCUCGACUCCUGGCC", 1817937833786320951),
            ## ("GGCCACCCUUGACCGUUUCCUCGUUCCGCUAGCCAACUUACGUGGUGGCAGAUUGGACGGUCAACUAUAUGUGGACAACUAUAAAUUAGUGGAGCGUAGCGACGGAUCUGUAUUGGUACCAGUGU+CCAAACCAGGCUGUCUCGACUCCUGGCC", 15606250151430488623),
            ## ("GGCCACCCUUGACCGUUUCCUCGUUCCGCUAGCCAACUUACGUGGUGGCAGAUUGGACGGUCAACUAUAUGUGGACAACUAUAAAUUAGUGGAGCGUAGCGACGGUCUGUUGGUACCAGUGU+CCAACAGGCUGUCCGCUCCUGGCC", 17872895797936292537),
            ## ("UAGGAAGAUCGAGUGGCGGCUUCAACCAGGUCGACCCUGUAAAGGGAAGCUUAAAAUCAAGGACCUACCG+UAUGCGACGGUCGCAUGACUUGGCGCUGAGCCAGUCCGCACUUUCGUAAUCGUAACUACGCACUAAGCGAAUACGUGGCGAGCCGGA", 2623438644795693868),
            ## ("GAAUACAGACAGAGUUAUGAAGGACGCUACGCAUUGUUAAAACGUCGUCUCGUUCAUUGGUCUCGGUUUAUUACGAAAGGAUGUCGGGUGCUUAAAGGGGGCUCACGUCUGACGCCACGCUCAUGGAGCACGCUGGGAUACAUAUUCCCAAUUGUUAGGGCCUGGAAGAC+CGGACCGACCAAGCAU", 15187308047988240186)
        ])

    def tearDown(self):
        shutil.rmtree(self.fp_workdir)  # remove temporary directory

    # tests unambiguoisness by checking if every dotBracket annotation appears exactly once
    def testUnambiguousness(self):
        os.system('basedir=`pwd` && cd "%s" && gapc -p "alg_cofold_dotBracket_unique*alg_count" $basedir/../../../cofold_nodangle.gap -I $basedir/../../../' % self.fp_workdir)
        os.system('basedir=`pwd` && cd "%s" && make -f out.mf CPPFLAGS_EXTRA="-I $basedir/../../../"' % self.fp_workdir)
        for seq in list(self.sequences.keys())[0:19]:
            cmd = ('%s/out "%s"' % (self.fp_workdir, seq))
            dot_count = [(subprocess.check_output(cmd, shell=True)).decode("utf-8")]
            print('  Testing with sequence "%s": got %i output lines.' % (seq, len(dot_count[0].split('\n'))), file=sys.stderr)
            for element in dot_count:
                dot_count_list = element.split('\n')
                for el in dot_count_list:
                    if not (el == 'Answer: ' or el == ''):
                        self.assertEqual(el.split(",")[1], ' 1 )')

    # tests if the predicted 2nd structures match a specific count 
    def testCount(self):
        os.system('basedir=`pwd` && cd "%s" && gapc -p "alg_count" $basedir/../../../cofold_nodangle.gap -I $basedir/../../../' % self.fp_workdir)
        os.system('basedir=`pwd` && cd "%s" && make -f out.mf CPPFLAGS_EXTRA="-I $basedir/../../../"' % self.fp_workdir)
        for seq, count in self.sequences.items():
            cmd = ('%s/out "%s"' % (self.fp_workdir, seq))
            candidates = (subprocess.check_output(cmd, shell=True)).decode("utf-8")
            print('  Testing with sequence "%s": got %s structures' % (seq, candidates.split('\n')[1]), file=sys.stderr)
            self.assertEqual(count[0], candidates.split('\n')[1])

    # tests mfe value
    def testMfe(self):
        os.system('basedir=`pwd` && cd "%s" && gapc -p "alg_cofold_mfe" $basedir/../../../cofold_nodangle.gap -I $basedir/../../../' % self.fp_workdir)
        os.system('basedir=`pwd` && cd "%s" && make -f out.mf CPPFLAGS_EXTRA="-I $basedir/../../../"' % self.fp_workdir)
        for seq, count in self.sequences.items():
            cmd = ('%s/out "%s"' % (self.fp_workdir, seq))
            mfe = (subprocess.check_output(cmd, shell=True)).decode("utf-8")
            print('  Testing with sequence "%s": got mfe of %s' % (seq, mfe.split('(')[1].split(',')[0]), file=sys.stderr)
            self.assertEqual(count[1], mfe.split('(')[1].split(',')[0])



if __name__ == "__main__":
    unittest.main()
