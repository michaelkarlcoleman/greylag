'''small run test cases'''


from __future__ import with_statement

import contextlib
import filecmp
import inspect
import subprocess

from nose.tools import *

from greylag import *

# program and initial args required to run from test directory
GREYLAG_PROGRAM = 'python2.5 ../greylag.py'

# temporary output file
GREYLAG_OUTPUT = 'tmp-greylag-out.xml'

# cd to the test directory, but restore when we're done
SAVE_CWD = os.getcwd()
def setup(self):
    if os.path.exists(GREYLAG_OUTPUT):
        os.remove(GREYLAG_OUTPUT)
    os.chdir('test')
def teardown(self):
    if os.path.exists(GREYLAG_OUTPUT):
        os.remove(GREYLAG_OUTPUT)
    os.chdir(SAVE_CWD)


def run_gl(args):
    "Run greylag in a subprocess, checking for error return."
    subprocess.check_call(("%s -q -o %s %s"
                           % (GREYLAG_PROGRAM, GREYLAG_OUTPUT, args)).split())

def run_combination(combination=None):
    """Run a greylag test, as specified by combination.

    If combination is not given, it defaults to the name of the calling
    function.

    For a combination 'greylag_params_0__test_2_test', greylag will be run
    with a parameter file 'greylag-params-0.xml' and spectrum file
    'test-2.ms2'.  (Additional spectrum files can be specified, separated by
    '__'.)  The results will be compared to 'greylag-params-0--test-2-ok.xml'.
    If this baseline file is not present and the environment variable
    NOSEUPDATE is set, a new baseline will be created.
    """
    if combination == None:
        # name of caller
        combination = inspect.getouterframes(inspect.currentframe())[1][3]
    assert combination.endswith('_test')
    combination = combination.rpartition('_test')[0]
    parts = combination.split('__')
    params, spectra = parts[0], parts[1:]
    assert len(spectra) >= 1
    params = params.replace('_', '-') + '.xml'
    spectra = [ sp.replace('_', '-') + '.ms2' for sp in spectra ]
    correct_fn = combination.replace('_', '-') + '-ok.xml'
    assert os.path.exists(correct_fn) or 'NOSEUPDATE' in os.environ, 'no baseline'
    run_gl(params + ' ' + ' '.join(spectra))
    if os.path.exists(correct_fn):
        assert filecmp.cmp(GREYLAG_OUTPUT, correct_fn), 'output differs'
    else:
        os.rename(GREYLAG_OUTPUT, correct_fn)
    # FIX: now do it in parts, and compare that result
    # also do it for --quirks?

# Could do xtandem and rough compare to quirks?  Is this worth doing, since we
# can expect at least small differences?  Probably we should have a separate
# framework for comparing and contrasting output from any two of { greylag,
# xtandem, sequest }



class modless_run_tests:
    def greylag_params_0__test_2_test(self):
        run_combination()
    def greylag_params_0__6323840_test(self):
        run_combination()
