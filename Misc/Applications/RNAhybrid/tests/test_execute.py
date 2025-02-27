import sys
import os
from glob import glob
from os.path import dirname, join
from io import StringIO
sys.path.append(dirname(__file__) + '/../')

from unittest import TestCase, main

from parse_gapc import *
from execute import complement, compose_call, execute_subprocess, cache_execute

def getFP(filepath):
    return join(dirname(__file__), filepath)

def _clean_cachefiles():
    # clean potential old cache files
    for fp_old in glob(join(os.environ["TMPDIR"], '*.cache_unittest')):
        os.remove(fp_old)

class TestExecute(TestCase):
    def setUp(self):
        _clean_cachefiles()

    def tearDown(self):
        _clean_cachefiles()

    def test_complement(self):
        exp = "GGGUUUCCC"
        obs = complement("CCCaaaGGG")
        self.assertEqual(exp, obs)

        exp = "AAAAAA"
        obs = complement("uuuttt")
        self.assertEqual(exp, obs)

        with self.assertRaises(TypeError):
            obs = complement("aaaya")

    def test_compose_call(self):
        exp = "./rnahybrid_eval_microstate -u 1 'CCCAAAGGG' '(((...)))'"
        obs = compose_call('eval', 'microstate', "CCCaaaGGG", None, inp_structure="(((...)))",
                           **{'binpath': './', 'program_name': 'rnahybrid', 'allow_lonelypairs': True})
        self.assertEqual(exp, obs)

        settings = {'binpath': './', 'program_name': 'rnahybrid'}
        exp = "./rnahybrid_khorshid_rnahybrid 'CCCCC' 'GGGGG'"
        obs = compose_call('khorshid', 'rnahybrid', "ccccc", "ggggg", **settings)
        self.assertEqual(exp, obs)

        exp = "./rnahybrid_mde 'UUCCC' 'GGGGG'"
        obs = compose_call('mde', '', 'tuCCC', complement('CCCCC'), **settings)
        self.assertEqual(exp, obs)

        with self.assertRaises(ValueError):
            compose_call('foo', '', 'CCCCC', complement('CCCCC'), **settings)

        with self.assertRaises(ValueError):
            compose_call('mde', 'foo', 'CCCCC', complement('CCCCC'), **settings)

    def test_execute_subprocess(self):
        exp = ['heinz', '']
        errFile = StringIO("")
        obs = execute_subprocess("echo 'heinz'", verbose=errFile)
        self.assertEqual(exp, obs)
        self.assertEqual("Binary call: echo 'heinz'\n", errFile.getvalue())

        with self.assertRaises(ValueError):
            obs = execute_subprocess("touch /proc", verbose=None)

    def test_cache_execute(self):
        cmd = "echo 'heinz'"
        errFile1 = StringIO("")
        # first execution, should create cache file
        obs1 = cache_execute(cmd, True, '.cache_unittest', verbose=errFile1)
        fp_cache = errFile1.getvalue().split()[-1].split("'")[1]
        self.assertEqual(['heinz', ''], obs1)
        self.assertIn('Wrote results into cache file', errFile1.getvalue())
        with open(fp_cache, 'r') as f:
            self.assertEqual('heinz\n', ''.join(f.readlines()))

        # second execution, should read from cache file. To test, we sneak a new value into the cache file
        with open(fp_cache, 'w') as f:
            f.write('cache is used\n')
        errFile2 = StringIO("")
        obs2 = cache_execute(cmd, True, '.cache_unittest', verbose=errFile2)
        self.assertEqual(['cache is used'], obs2)
        self.assertIn('Read cached result from file', errFile2.getvalue())

        # remove cache files
        _clean_cachefiles()

        # execution without using a cache file
        obs3 = cache_execute(cmd, False, '.cache_unittest', verbose=None)
        self.assertEqual(['heinz', ''], obs3)


if __name__ == '__main__':
    main()
