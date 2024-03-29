'''small run test cases'''


from __future__ import with_statement

import contextlib
import filecmp
import inspect
import subprocess

from nose.tools import *

from greylag_chase import *


# greylag-solo parallel processes
CPUS = 4

# program locations
GREYLAG_CHASE = '../greylag_chase.py'
GREYLAG_SOLO = '../greylag_solo.py'
GREYLAG_SQT = '../greylag_sqt.py'
GREYLAG_INDEX_SPECTRA = '../greylag_index_spectra.py'

JOB_ID = 'tmp'
DEFAULT_GLW = 'chase_%s_0.0-1.0.glw' % JOB_ID

# cd to the test directory, but restore when we're done
SAVE_CWD = os.getcwd()
def setup(self):
    os.chdir('test')
def teardown(self):
    os.chdir(SAVE_CWD)


def run_gl(conf, spectra):
    "Run greylag-chase, checking for error return."
    try:
        os.remove(DEFAULT_GLW)
    except:
        pass
    subprocess.check_call([GREYLAG_CHASE, "-q", "-w", "0:1",
                           "--job-id="+JOB_ID, conf] + spectra)

def run_gl_solo(conf, spectra):
    "Run greylag-solo, checking for error return."
    subprocess.check_call([GREYLAG_SOLO, str(CPUS), "-q", conf] + spectra)

def run_gl_sqt(arg):
    "Run greylag-sqt, checking for error return."
    subprocess.check_call([GREYLAG_SQT, arg])

def run_gl_index_spectra(args):
    "Run greylag-index-spectra, checking for error return."
    subprocess.check_call([GREYLAG_INDEX_SPECTRA] + args)


# This is currently kind of bogus because ../greylag_solo.py will end up
# running the installed versions of greylag-chase, etc.  How to fix this?

def run_combination(combination=None):
    """Run a greylag test, as specified by combination.

    If combination is not given, it defaults to the name of the calling
    function.

    For a combination 'greylag_params_0__test_2_test', greylag will be run
    with a parameter file 'greylag-params-0.conf' and spectrum file
    'test-2.ms2'.  (Additional spectrum files can be specified, separated by
    '__'.)  The results will be compared to 'greylag-params-0--test-2-ok.sqt'.
    If this baseline file is not present and the environment variable
    NOSEUPDATE is set, a new baseline will be created.
    """

    if combination == None:
        # name of caller
        combination = inspect.getouterframes(inspect.currentframe())[1][3]
    assert combination.endswith('_test')
    combination = id = combination.rpartition('_test')[0]
    run = run_gl
    if combination.endswith('_solo'):
        run = run_gl_solo
        combination = combination.rpartition('_solo')[0]
    parts = combination.split('__')
    params, spectra = parts[0], parts[1:]
    assert len(spectra) >= 1
    params = params.replace('_', '-') + '.conf'
    ms2 = [ sp.replace('_', '-') + '.ms2' for sp in spectra ]
    sqt = [ sp.replace('_', '-') + '.sqt' for sp in spectra ]
    sqt_ok = [ sp.replace('_', '-') + '-ok.sqt' for sp in spectra ]
    if not 'NOSEUPDATE' in os.environ:
        for sok in sqt_ok:
            assert os.path.exists(sok), "baseline '%s' missing" % sok

    run_gl_index_spectra(ms2)
    run(params, ms2)
    if run == run_gl:
        run_gl_sqt(DEFAULT_GLW)
        os.remove(DEFAULT_GLW)

    for s, sok in zip(sqt, sqt_ok):
        assert os.path.exists(s), "sqt output '%s' missing" % s
        if os.path.exists(sok):
            if not filecmp.cmp(s, sok):
                os.rename(s, "%s_%s-bad.sqt" % (os.path.splitext(s)[0], id))
                assert False, "sqt output '%s' differs" % s
            os.remove(s)
        else:
            os.rename(s, sok)


# FIX
# class modless_run_tests:
#     def greylag_params_0__test_2_test(self):
#         run_combination()
#     def greylag_params_0__test_2_solo_test(self):
#         run_combination()
#     def greylag_params_0__6323840_test(self):
#         run_combination()
#     def greylag_params_0__6323840_solo_test(self):
#         run_combination()

#     # these take forever, unfortunately
#     # def greylag_params_yeast_0__test_2_test(self):
#     #     run_combination()
#     # def greylag_params_yeast_0__test_2_solo_test(self):
#     #     run_combination()


# class methox_run_tests:
#     def methox__ScTAP_Ceg1_YPD_p1_Ti_108_007713_test(self):
#         run_combination()
#     def methox__ScTAP_Ceg1_YPD_p1_Ti_108_007713_solo_test(self):
#         run_combination()

