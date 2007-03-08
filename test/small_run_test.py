'''small run test cases'''


from __future__ import with_statement

import contextlib
import filecmp
import subprocess

from nose.tools import *

from greylag import *

# program and initial args required to run from test directory
GREYLAG_PROGRAM = 'python2.5 ../greylag.py'

# temporary greylag
GREYLAG_OUTPUT = 'tmp-greylag-out.xml'


# cd to the test directory, but restore when we're done
SAVE_CWD = os.getcwd()
def setup(self):
    os.chdir('test')
def teardown(self):
    os.chdir(SAVE_CWD)


def run_gl(args):
    "Run greylag in a subprocess, checking for error return."
    subprocess.check_call(("%s -o %s %s"
                           % (GREYLAG_PROGRAM, GREYLAG_OUTPUT, args)).split())


class modless_run_tests:
    def greylag_params_0_test(self):
        run_gl('-q greylag-params-0.xml test-2.ms2')
        assert filecmp.cmp(GREYLAG_OUTPUT, 'greylag-params-0-ok.xml')
