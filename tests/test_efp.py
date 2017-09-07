#import pylibefp
import molparse
from utils import *
from addons import *


_efpefp1_str = """
efp C6H6 -0.30448173 -2.24210052 -0.29383131 -0.642499 1.534222 -0.568147
--
efp C6H6 -0.60075437  1.36443336  0.78647823  3.137879 1.557344 -2.568550
"""


_efpefp1_rdict = {'libefp': {'full_fragments': 
 [{'coordinates_hint': [-0.5753870821672306,
                       -4.23695594520049,
                       -0.5552607051670226,
                       -0.642499,
                       1.534222,
                       -0.568147],
  'efp_type': 'xyzabc',
  'fragment_file': 'C6H6'},
 {'coordinates_hint': [-1.1352612324342508,
                       2.578405376972965,
                       1.4862284641766452,
                       3.137879,
                       1.557344,
                       -2.56855],
  'efp_type': 'xyzabc',
  'fragment_file': 'C6H6'}]},
 'molecule': {
 'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [],
 'fragment_multiplicities': [],
 'fragment_types': [],
 'fragments': [],
 'full_atoms': [],
 'input_units_to_au': 1.8897261328856432,
 'name': 'default',
 'units': 'Angstrom'}}


def test_efpefp1():
    """qchem-efp-sp"""

    pdict = molparse.init_psi4_molecule_from_any_string(_efpefp1_str)
    #print('RDICT', _efpefp1_rdict, '\n\n', 'WDICT', pdict, '\n')
    assert(compare_dicts(_efpefp1_rdict, pdict, 6, sys._getframe().f_code.co_name + ' str --> dict'))


@using_qcdb
def test_efpefp1_q():
    import qcdb
    pdict = molparse.init_psi4_molecule_from_any_string(_efpefp1_str)

    rnre = 0.000

    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
    assert(compare_values(rnre, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ' dict --> qcdbobj'))

    qdict = qmol.to_dict()
    print('RDICT', _efpefp1_rdict, '\n\n', 'QDICT', qdict, '\n')
    assert(compare_dicts(_efpefp1_rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ' qcdbobj --> dict'))
    # todo need another for _efpefp1_rdict['libefp']


_libefp1_str = """
efp h2o 0.0 0.0 0.0 1.0, 2.0, 3.0
--
efp nh3 5.0 0.0 0.0 5.0, 2.0, 8.0
"""


_libefp1_rdict = {'libefp': {'full_fragments':
 [{'coordinates_hint': [0.0,
                       0.0,
                       0.0,
                       1.0,
                       2.0,
                       3.0],
  'efp_type': 'xyzabc',
  'fragment_file': 'h2o'},
 {'coordinates_hint': [9.448630664428215,
                       0.0,
                       0.0,
                       5.0,
                       2.0,
                       8.0],
  'efp_type': 'xyzabc',
  'fragment_file': 'nh3'}]},
 'molecule': {
 'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [],
 'fragment_charges': [],
 'fragment_multiplicities': [],
 'fragment_types': [],
 'fragments': [],
 'full_atoms': [],
 'input_units_to_au': 1.8897261328856432,
 'name': 'default',
 'units': 'Angstrom'}}

_libefp1_vdict = {'libefp': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0, 0],
 'fragment_multiplicities': [1, 1],
 'fragment_types': ['Real', 'Real'],
 'fragments': [[0, 3], [3, 7]],
 'full_atoms': [{'Z': 8.0,
                 'charge': 8.0,
                 'ghosted': False,
                 'label': 'A01O1',
                 'mass': 15.99491461956,
                 'at_type': 'efpxyz',
                 'symbol': 'O',
                 'x': -0.050373610528,
                 'y': 0.012369110219,
                 'z': -0.107222070225},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'A02H2',
                 'mass': 1.00782503207,
                 'at_type': 'efpxyz',
                 'symbol': 'H',
                 'x': 1.090231686178,
                 'y': 1.131828547911,
                 'z': 0.668335425557},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'A03H3',
                 'mass': 1.00782503207,
                 'at_type': 'efpxyz',
                 'symbol': 'H',
                 'x': -0.290766137580,
                 'y': -1.328135252620,
                 'z': 1.033356200180},
                {'Z': 7.0,
                 'charge': 7.0,
                 'ghosted': False,
                 'label': 'A01N1',
                 'mass': 14.00307400478,
                 'at_type': 'efpxyz',
                 'symbol': 'N',
                 'x': 9.552734528021,
                 'y': 0.030794165669,
                 'z': 0.049682979512},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'A02H2',
                 'mass': 1.00782503207,
                 'at_type': 'efpxyz',
                 'symbol': 'H',
                 'x': 9.691240950238,
                 'y': -1.638423923013,
                 'z': -0.821527157536},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'A03H3',
                 'mass': 1.00782503207,
                 'at_type': 'efpxyz',
                 'symbol': 'H',
                 'x': 9.011097119400,
                 'y': 1.313072001329,
                 'z': -1.225806192962},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'A04H4',
                 'mass': 1.00782503207,
                 'at_type': 'efpxyz',
                 'symbol': 'H',
                 'x': 8.197161749192,
                 'y': -0.102512891040,
                 'z': 1.357020635862}],
 'input_units_to_au': 1.0,
 'units': 'Bohr'}}
# 'input_units_to_au': 1.8897261328856432,
# 'name': 'default',
# 'units': 'Angstrom'}}


def test_libefp1_b():

    pdict = molparse.init_psi4_molecule_from_any_string(_libefp1_str)
    print('RDICT', _libefp1_rdict, '\n\n', 'WDICT', pdict, '\n')
    assert(compare_dicts(_libefp1_rdict, pdict, 6, sys._getframe().f_code.co_name + ': [1] str --> dict'))


@using_qcdb
def test_libefp1_q():
    import qcdb
    pdict = molparse.init_psi4_molecule_from_any_string(_libefp1_str)

    rnreqm = 0.000

    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
    assert(compare_values(rnreqm, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': [2] dict --> qcdbobj'))

    qdict = qmol.to_dict()
    #print('RDICT', _libefp1_rdict['molecule'], '\n\n', 'QDICT', qdict, '\n')
    assert(compare_dicts(_libefp1_rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ': [3] qcdbobj --> dict'))


@using_pylibefp
def test_libefp1_e():
    import pylibefp

    pdict = molparse.init_psi4_molecule_from_any_string(_libefp1_str)

    rnreefp = 32.128

    emol = pylibefp.from_dict(pdict['libefp'])
    assert(compare_values(rnreefp, emol.nuclear_repulsion_energy(), 3, sys._getframe().f_code.co_name + ': [4] dict --> efpobj'))

    #edict = emol.to_dict()
    edict = emol.to_viz_dict()
    print('\nRDICT')
    pprint.pprint(_libefp1_vdict['libefp'])
    print('\nEDICT')
    pprint.pprint(edict)
#    print('RDICT', _libefp1_vdict['libefp'], '\n\n', 'EDICT', edict, '\n')
    assert(compare_dicts(_libefp1_vdict['libefp'], edict, 4, sys._getframe().f_code.co_name + ': [5] efpobj --> dict'))
    # _libefp1_rdict['libefp']

    

_qmefp1_str = """
# QM fragment
0 1
O    0         0.0   0.118720
H   -0.753299, 0.0, -0.474880
H    0.753299, 0.0, -0.474880
# EFP as EFP fragments
--
efp water
    -2.13972    1.28964   -0.96418
    -2.66865    0.51034   -1.14473
    -1.33300    0.93113   -0.58956
--
efp ammonia
     0.98792  ,  1.87681    2.85174
1.68798,,1.18856    3.09517
     1.45873    2.55904    2.27226
--
efp ammonia
    -4.12794   -0.92466   -1.28394
    -4.69278   -1.09557   -2.10539
    -3.59191   -1.76923   -1.13470
"""


_qmefp1_rdict = {'libefp': {'full_fragments':
 [{'coordinates_hint': [-4.043484801,   2.437066410,  -1.822036143,
                        -5.043016189,   0.964404979,  -2.163225699,
                        -2.519001456,   1.759578808,  -1.114105301],
  'efp_type': 'points',
  'fragment_file': 'water'},
  {'coordinates_hint': [ 1.866898241,   3.546656903,   5.389007602,
                         3.189822098,   2.246050749,   5.849024393,
                         2.756589676,   4.835894818,   4.293945633],
  'efp_type': 'points',
  'fragment_file': 'ammonia'},
  {'coordinates_hint': [-7.800676093,  -1.747354166,  -2.426294971,
                        -8.868068505,  -2.070327109,  -3.978609780,
                        -6.787723797,  -3.343360142,  -2.144269196],
  'efp_type': 'points',
  'fragment_file': 'ammonia'}]},
 'molecule': {
 'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0],
 'fragment_multiplicities': [1],
 'fragment_types': ['Real'],
 'fragments': [[0, 3]],
 'full_atoms': [{'Z': 8.0,
                 'charge': 8.0,
                 'ghosted': False,
                 'label': 'O',
                 'mass': 15.99491461956,
                 'qm_type': 'qmcart',
                 'symbol': 'O',
                 'x':  0.0,         'y':  0.0,         'z':  0.224348286},
                 ##'x':  0,           'y':  0,           'z':  0.224348286},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': -1.423528806, 'y':  0.0,         'z': -0.897393146},
                 ##'x': -1.423528806, 'y':  0,           'z': -0.897393146},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':  1.423528806, 'y':  0.0,         'z': -0.897393146}],
                 ##'x':  1.423528806, 'y':  0,           'z': -0.897393146}],
 #'input_units_to_au': 1.8897261328856432,
 'name': 'default',
 'system_charge': 0,
 'system_multiplicity': 1,
 #'units': 'Angstrom'}}
 'input_units_to_au': 1.0,
 'units': 'Bohr'}}


_qmefp1_vdict = {'libefp': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0, 0, 0, 0],
 'fragment_multiplicities': [1, 1, 1, 1],
 'fragment_types': ['Real', 'Real', 'Real', 'Real'],
 'fragments': [[0, 3], [3, 6], [6, 10], [10, 14]],
 'full_atoms': [{'charge': 8.0, 'ghosted': False, 'label': '0', 'mass': 15.99491461956, 'at_type': 'qmxyz', 'symbol': 'O', 'Z': 8,   'x':  0,           'y':  0,           'z':  0.224348286},
                {'charge': 1.0, 'ghosted': False, 'label': 'H_Aa', 'mass': 1.00782503207, 'at_type': 'qmxyz', 'symbol': 'H', 'Z': 1,   'x': -1.423528806, 'y':  0,           'z': -0.897393146},
                {'charge': 1.0, 'ghosted': False, 'label': 'H_Bb', 'mass': 1.00782503207, 'at_type': 'qmxyz', 'symbol': 'H', 'Z': 1,   'x':  1.423528806, 'y':  0,           'z': -0.897393146},
                {'charge': 8.0, 'ghosted': False, 'label': 'A01O1', 'mass': 15.99491461956, 'at_type': 'efpxyz', 'symbol': 'O', 'Z': 8,   'x': -4.043484801, 'y':  2.437066410, 'z': -1.822036143},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x': -5.043016189, 'y':  0.964404979, 'z': -2.163225699},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x': -2.519001456, 'y':  1.759578808, 'z': -1.114105301},
                {'charge': 7.0, 'ghosted': False, 'label': 'A01N1', 'mass': 14.00307400478, 'at_type': 'efpxyz', 'symbol': 'N', 'Z': 7,   'x':  1.866898241, 'y':  3.546656903, 'z':  5.389007602},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x':  3.189822098, 'y':  2.246050749, 'z':  5.849024393},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x':  2.756589676, 'y':  4.835894818, 'z':  4.293945633},
                {'charge': 1.0, 'ghosted': False, 'label': 'A04H4', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x':  1.431508879, 'y':  4.450724844, 'z':  7.015643157},
                {'charge': 7.0, 'ghosted': False, 'label': 'A01N1', 'mass': 14.00307400478, 'at_type': 'efpxyz', 'symbol': 'N', 'Z': 7,   'x': -7.800676093, 'y': -1.747354166, 'z': -2.426294971},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x': -8.868068505, 'y': -2.070327109, 'z': -3.978609780},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x': -6.787723797, 'y': -3.343360142, 'z': -2.144269196},
                {'charge': 1.0, 'ghosted': False, 'label': 'A04H4', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1,   'x': -9.023788093, 'y': -1.636376864, 'z': -0.961874594}],
 'input_units_to_au': 1.0,
 'units': 'Bohr'}}


def test_qmefp1_b():
    """qchem-qmefp-puream-sp"""

    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp1_str)
    print('RDICT', _qmefp1_rdict, '\n\n', 'PDICT', pdict, '\n')
    assert(compare_dicts(_qmefp1_rdict, pdict, 3, sys._getframe().f_code.co_name + ': [1] str --> dict'))


#@using_qcdb
#def test_qmefp1_q():
#    import qcdb
#    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp1_str)
#
#    rnreqm = 0.000
#
#    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
#    assert(compare_values(rnreqm, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': [2] dict --> qcdbobj'))
#
#    qdict = qmol.to_dict()
#    #print('RDICT', rdict['molecule'], '\n\n', 'QDICT', qdict, '\n')
#    assert(compare_dicts(rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ': [3] qcdbobj --> dict'))
#
#
#@using_pylibefp
#def test_qmefp1_e():
#    import pylibefp
#    rnreefp = 32.128
#
#    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp1_str)
#
#    emol = pylibefp.from_dict(pdict['libefp'])
#    assert(compare_values(rnreefp, emol.nuclear_repulsion_energy(), 3, sys._getframe().f_code.co_name + ': [4] dict --> efpobj'))
#
#    edict = emol.to_dict()
#    #print('RDICT', rdict['libefp'], '\n\n', 'EDICT', edict, '\n')
#    assert(compare_dicts(_qmefp1_vdict['libefp'], edict, 6, sys._getframe().f_code.co_name + ': [5] efpobj --> dict'))
#    # rdict['libefp']


if __name__ == "__main__":
    test_efpefp1()
    #test_libefp1()
    #test_qmefp1()

    test_efpefp1()
    test_efpefp1_q()
    test_libefp1_b()
    test_libefp1_q()
    test_libefp1_e()
    test_qmefp1_b()
    test_qmefp1_q()
    test_qmefp1_e()

