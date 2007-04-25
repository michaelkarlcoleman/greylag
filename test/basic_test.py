'''simple, quick unit test cases'''


from __future__ import with_statement

import contextlib

from nose.tools import *

from greylag_grind import *


def same_atoms_test():
    assert (set(MONOISOTOPIC_ATOMIC_MASS.keys())
            == set(AVERAGE_ATOMIC_MASS.keys()))


class formula_mass_test:
    def summation_test(self):
        a = formula_mass('C5H9ONS')
        b = formula_mass('C5H9   ')
        c = formula_mass('    ONS')
        assert_almost_equal(a, b+c, places=10)


class read_fasta_files_test:
    def simple_test_0(self):
        found = list(read_fasta_files(['test/scary.fa']))
        wanted = [('ugly test fasta file', '', 'test/scary.fa'),
                  ('', 'ABC', 'test/scary.fa'),
                  ('1', 'STALEMATING', 'test/scary.fa'),
                  (' foo', 'AARDVARKS AARDVARKSAARDVARKS', 'test/scary.fa'),
                  ('gi|1234| bleah', 'SISTERLIKE', 'test/scary.fa'),
                  ('gi|1234| bleah', 'SISTER*LIKE', 'test/scary.fa'),
                  ('control-A here: [\x01]', '& HERE: [\x01]', 'test/scary.fa'),
                  ('no final newline after F', 'ARF', 'test/scary.fa')]
        assert found == wanted


def read_taxonomy_test():
    found = sorted(read_taxonomy('test/test-taxonomy.xml').items())
    wanted = [('label', ['labeled.fa']),
              ('multiple', ['one.fa', 'two.fa', 'three.fa']),
              ('simple.fa', ['simple.fa'])]
    assert found == wanted


class read_xml_parameters_test:
    def empty_test(self):
        assert [] == read_xml_parameters('test/empty-params.xml').items()


class mass_regime_part_test:
    @raises(ValueError)
    def bad_regime_test(self):
        mass_regime_part('NOT_MONO_OR_AVG')

    @raises(ValueError)
    def bad_syntax_0_test(self):
        mass_regime_part('MONO(N15@100%')

    @raises(ValueError)
    def bad_syntax_1_test(self):
        mass_regime_part('MONO(N15)')

    @raises(ValueError)
    def bad_syntax_2_test(self):
        mass_regime_part('MONO(N15@100)')

    @raises(ValueError)
    def N15_test(self):
        mass_regime_part('MONO(C13@99%)')

    @raises(ValueError)
    def multiple_isotope_test(self):
        mass_regime_part('MONO(N15@80%,C13@99%)')

    @raises(ValueError)
    def bad_prevalence_0_test(self):
        mass_regime_part('MONO(N15@100.1%)')

    @raises(ValueError)
    def bad_prevalence_1_test(self):
        mass_regime_part('MONO(N15@-1.1%)')


class mass_regime_list_test:
    @raises(ValueError)
    def bad_frag_test(self):
        mass_regime_list('AVG/AVG')


class parse_mod_term_test:
    @raises(ValueError)
    def bad_spec_0_test(self):
        parse_mod_term('-C2H3ON!!@C')

    @raises(ValueError)
    def bad_spec_1_test(self):
        parse_mod_term('-2C2H3ON!@C')

    @raises(ValueError)
    def bad_spec_2_test(self):
        parse_mod_term('-2Z@C')

    @raises(ValueError)
    def bad_spec_4_test(self):
        parse_mod_term('--C2H3ON@C')

    @raises(ValueError)
    def bad_spec_5_test(self):
        parse_mod_term('-C2H3ON@C foo bar')

    @raises(ValueError)
    def bad_spec_6_test(self):
        parse_mod_term('-C2H3ON@C@G')

    @raises(ValueError)
    def bad_spec_7_test(self):
        parse_mod_term('-C2H3ON@@C')

    @raises(ValueError)
    def bad_spec_8_test(self):
        parse_mod_term('++C2H3ON@C')

    @raises(ValueError)
    def bad_spec_9_test(self):
        parse_mod_term('+C2H3ON@C "foo"')

    @raises(ValueError)
    def bad_residue_0_test(self):
        parse_mod_term('20@Z')

    @raises(ValueError)
    def bad_residue_1_test(self):
        parse_mod_term('20@[C', is_potential=True)

    @raises(ValueError)
    def bad_residue_2_test(self):
        parse_mod_term('20@]C', is_potential=True)

    @raises(ValueError)
    def bad_delta_0_test(self):
        parse_mod_term('0@C', is_potential=True)

    @raises(ValueError)
    def bad_delta_1_test(self):
        parse_mod_term('0.000000001@C', is_potential=True)

    @raises(ValueError)
    def duplicate_residue_0_test(self):
        parse_mod_term('120.0@CC')

    @raises(ValueError)
    def duplicate_residue_1_test(self):
        parse_mod_term('120.0@CC', is_potential=True)

    @raises(ValueError)
    def duplicate_residue_2_test(self):
        parse_mod_term('120.0@STY')


class fixed_mod_list_test:
    @raises(ValueError)
    def multiple_residue_test(self):
        fixed_mod_list('14@ST')


class parse_mod_basic_expression_test:
    @raises(ValueError)
    def missing_paren_0_test(self):
        parse_mod_basic_expression('(80@STY;12@C')

    @raises(ValueError)
    def missing_paren_1_test(self):
        parse_mod_basic_expression('(80@STY')

    @raises(ValueError)
    def missing_paren_2_test(self):
        parse_mod_basic_expression('(((80@STY))')




class potential_mod_list_test:
    @raises(ValueError)
    def extra_paren_0_test(self):
        potential_mod_list('(((80@STY;12@C))))')

    @raises(ValueError)
    def extra_paren_0_test(self):
        potential_mod_list('(((80@STY;12@C)))(')


def XML_PARAMETER_INFO_test():
    for k, v in XML_PARAMETER_INFO.iteritems():
        kv = "%s: %s" % (k, v)
        assert isinstance(k, str)
        assert isinstance(v, tuple)
        assert 2 <= len(v) <= 3
        if v[0] == bool:
            assert v[1] == None or v[1] in ('yes', 'no'), kv
        elif v[0] in (int, float, str):
            assert v[1] == None or isinstance(v[0](v[1]), v[0]), kv
        elif isinstance(v[0], tuple):
            assert v[1] == None or v[1] in v[0], kv
        elif v[0] in (mass_regime_list, fixed_mod_list, potential_mod_list):
            assert v[1] == None or isinstance(v[1], str), kv
        else:
            assert False, kv


class zopen_test:
    def read_check(self, filename):
        with contextlib.closing(zopen(filename)) as f:
            data = f.read()
        assert data == 'data\n'

    def plain_test(self):
        self.read_check('test/data')
    def gz_test(self):
        self.read_check('test/data.gz')
    def bz2_test(self):
        self.read_check('test/data.bz2')


class swig_sanity_test:
    def vector_test(self):
        x = cgreylag.mass_regime_parameters()
        x.fixed_residue_mass.resize(100)
        assert len(x.fixed_residue_mass) == 100
        x.fixed_residue_mass.resize(200, 1.1)
        assert len(x.fixed_residue_mass) == 200
        for i in range(0, 100):
            assert x.fixed_residue_mass[i] == 0.0
        for i in range(100, 200):
            assert x.fixed_residue_mass[i] == 1.1


def reset_spectrum_ids():
    "reset spectrum serial numbers, so that tests are repeatable"
    cgreylag.cvar.spectrum_next_id = 0
    cgreylag.cvar.spectrum_next_physical_id = 0


class external_type_test:
    def peak_repr_test(self):
        x = cgreylag.peak()
        assert repr(x) == '<peak mz=0.0000 intensity=0.0000>'
        x = cgreylag.peak(1234.56, 112233)
        assert repr(x) == '<peak mz=1234.5600 intensity=112233.0000>'

    def spectrum_repr_test(self):
        reset_spectrum_ids()
        x = cgreylag.spectrum()
        assert repr(x) == "<spectrum #0 (phys #-1) '' 0.0000/+0 [0 peaks, maxI=-1.000000, sumI=-1.000000]>"
        x = cgreylag.spectrum(1234.12, 2)
        assert repr(x) == "<spectrum #1 (phys #-1) '' 1234.1200/+2 [0 peaks, maxI=-1.000000, sumI=-1.000000]>"

    def set_peaks_from_matrix_test(self):
        reset_spectrum_ids()
        x = cgreylag.spectrum()
        x.set_peaks_from_matrix([(1.0,2.0), (3.5,4.5), (5.0,6.0)])
        assert repr(x) == "<spectrum #0 (phys #-1) '' 0.0000/+0 [3 peaks, maxI=-1.000000, sumI=-1.000000]>"



class spectrum_test:
    def simple_read_test(self):
        reset_spectrum_ids()
        FILEID = 8
        with open('test/simple.ms2') as f:
            x = cgreylag.spectrum.read_spectra_from_ms2(f, FILEID)
        assert len(x) == 6
        wanted = ("<spectrum #0 (phys #0) '0002.0002.2' 2096.1000/+2 [356 peaks, maxI=2475458.000000, sumI=19305560.000000]>",
                  "<spectrum #1 (phys #0) '0002.0002.3' 3143.6600/+3 [356 peaks, maxI=2475458.000000, sumI=19305560.000000]>",
                  "<spectrum #2 (phys #1) '0003.0003.2' 2307.4200/+2 [160 peaks, maxI=13492.000000, sumI=535952.000000]>",
                  "<spectrum #3 (phys #1) '0003.0003.3' 3460.6000/+3 [160 peaks, maxI=13492.000000, sumI=535952.000000]>",
                  "<spectrum #4 (phys #2) '0006.0006.2' 2057.9100/+2 [410 peaks, maxI=689515.000000, sumI=4639641.000000]>",
                  "<spectrum #5 (phys #2) '0006.0006.3' 3086.2400/+3 [410 peaks, maxI=689515.000000, sumI=4639641.000000]>")
        for n, sp in enumerate(x):
            assert repr(sp) == wanted[n], "checking %s" % n
            assert sp.file_id == FILEID

    @raises(RuntimeError)
    def ms2_instead_of_ms2_plus_test(self):
        reset_spectrum_ids()
        with open('test/simple.ms2') as f:
            x = cgreylag.spectrum.read_spectra_from_ms2(f, -1)


class read_spectra_test:
    F = 'test/junk.ms2'
    def setup(self):
        self.create_ms2([''], False)    # just checking that we can

    def create_ms2(self, contents, final_newline=True):
        assert not isinstance(contents, str)
        with open(self.F, 'w') as f:
            f.write('\n'.join(contents))
            if final_newline:
                f.write('\n')

    def read_ms2(self, contents, final_newline=True):
        reset_spectrum_ids()
        self.create_ms2(contents, final_newline)
        with open(self.F) as f:
            return cgreylag.spectrum.read_spectra_from_ms2(f, 1)

    def test_empty(self):
        assert len(self.read_ms2([''], False)) == 0

    def test_wordy_name(self):
        lines = [':blah blah blah', '1234.5 1', '123 456']
        wanted = "(<spectrum #0 (phys #0) 'blah blah blah' 1234.5000/+1 [1 peaks, maxI=456.000000, sumI=456.000000]>,)"
        assert wanted == repr(self.read_ms2(lines))

    @raises(RuntimeError)
    def test_no_peaks_at_eof(self):
        self.read_ms2([':0002.0002.1', '1234.5 1'])

    @raises(RuntimeError)
    def test_missing_mass_charge(self):
        self.read_ms2([':'])

    @raises(RuntimeError)
    def test_missing_additional_mass_charge(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', ':0002.0002.3'])

    @raises(RuntimeError)
    def test_extra_following_mass_charge(self):
        self.read_ms2([':0002.0002.2', '1234.5 2 3', '123 456 789'])

    @raises(RuntimeError)
    def test_extra_following_peak(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', '123 456 789'])

    @raises(RuntimeError)
    def test_missing_colon_line(self):
        self.read_ms2(['1234.5 2', '123 456 789'])

    @raises(RuntimeError)
    def test_bad_number_0(self):
        self.read_ms2([':0002.0002.2', 'a1234.5 2', '123 456'])

    @raises(RuntimeError)
    def test_bad_number_1(self):
        self.read_ms2([':0002.0002.2', '1234.5 a2', '123 456'])

    @raises(RuntimeError)
    def test_bad_number_2(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', 'a123 456'])

    @raises(RuntimeError)
    def test_bad_number_3(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', '123 a456'])

    @raises(RuntimeError)
    def test_nonpositive_number_0(self):
        self.read_ms2([':0002.0002.2', '-1234.5 2', '123 456'])

    @raises(RuntimeError)
    def test_nonpositive_number_1(self):
        self.read_ms2([':0002.0002.2', '1234.5 -2', '123 456'])

    @raises(RuntimeError)
    def test_nonpositive_number_2(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', '-123 456'])

    @raises(RuntimeError)
    def test_nonpositive_number_3(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', '123 -456'])

    @raises(RuntimeError)
    def test_zero_number_0(self):
        self.read_ms2([':0002.0002.2', '0.0 2', '123 456'])

    @raises(RuntimeError)
    def test_zero_number_1(self):
        self.read_ms2([':0002.0002.2', '1234.5 0', '123 456'])

    @raises(RuntimeError)
    def test_zero_number_2(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', '0.0 456'])

    @raises(RuntimeError)
    def test_bad_blank_0(self):
        self.read_ms2([''])

    @raises(RuntimeError)
    def test_bad_blank_1(self):
        self.read_ms2([':0002.0002.2', '', '1234.5 2', '123 456'])

    @raises(RuntimeError)
    def test_bad_blank_2(self):
        self.read_ms2([':0002.0002.2', '1234.5 2', '', '123 456'])

    def test_read_masses_0(self):
        ms = cgreylag.spectrum.read_ms2_spectrum_masses(())
        assert ms == ()

    def test_read_masses_1(self):
        self.create_ms2([':', '10234.5 2', '123 456',
                         ':', '30234.5 2', '123 456',
                         ':', '20234.5 2', '123 456',])
        with open(self.F) as f:
            fn = f.fileno()
            ms = cgreylag.spectrum.read_ms2_spectrum_masses((fn,))
        assert ms == (10234.5, 20234.5, 30234.5)

    def test_read_masses_2(self):
        self.create_ms2([':', '102.5 2', '123 456',
                         ':', '302.5 2', '123 456',
                         ':', '202.5 2', '123 456',])
        with open(self.F) as f1:
            fn1 = f1.fileno()
            with open(self.F) as f2:
                fn2 = f2.fileno()
                ms = cgreylag.spectrum.read_ms2_spectrum_masses((fn1, fn2))
                assert ms == (102.5, 102.5, 202.5, 202.5, 302.5, 302.5)

    def teardown(self):
        os.remove(self.F)


# FIX: more cgreylag testing to do here
