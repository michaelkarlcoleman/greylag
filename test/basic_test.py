'''simple, quick unit test cases'''


from __future__ import with_statement

import contextlib

from nose.tools import *

from greylag import *


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


class test_zopen:
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
