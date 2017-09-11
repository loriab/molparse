import molparse
from utils import *
from addons import *


def filter_radian_range(dicary):
    """Adjust any `dicary` xyzabc efp fragment hints by 2pi into the
    (-pi, pi] range.

    """
    try:
        asdf = dicary['libefp']['full_fragments']
    except KeyError:
        return dicary

    def radrge(elem):
        if elem > math.pi:
            return elem - 2 * math.pi
        elif elem <= -math.pi:
            return elem + 2 * math.pi
        else:
            return elem

    for idx, itm in enumerate(asdf):
        if itm['efp_type'] == 'xyzabc':
            hint = itm['coordinates_hint']
            hint[3] = radrge(hint[3])
            hint[4] = radrge(hint[4])
            hint[5] = radrge(hint[5])
            asdf[idx]['coordinates_hint'] = hint
    dicary['libefp']['full_fragments'] = asdf
    return dicary


_efpefp1_nre = 635.898999423     #      EFP(Bz-Bz)
_libefp1_nre = 32.1285883797     #      EFP(NH3-H2O)
_qmefp1_nre_tot = 124.46435096   #      EFP(H2O-NH3-NH3)-QM(H2O)
_qmefp1_nre_qm = 9.17938788302   #  QM@ EFP(H2O-NH3-NH3)-QM(H2O)
_qmefp1_nre_efp = 68.2120352668  # EFP@ EFP(H2O-NH3-NH3)-QM(H2O)
_qmefp2_nre_tot = 370.106323977  #      EFP(H2O-CL2)-QM(H2O)
_qmefp2_nre_qm = 8.35355339059   #  QM@ EFP(H2O-CL2)-QM(H2O)
_qmefp2_nre_efp = 115.631159467  # EFP@ EFP(H2O-CL2)-QM(H2O)


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
  'fragment_file': 'c6h6'},
 {'coordinates_hint': [-1.1352612324342508,
                       2.578405376972965,
                       1.4862284641766452,
                       3.137879,
                       1.557344,
                       -2.56855],
  'efp_type': 'xyzabc',
  'fragment_file': 'c6h6'}],
 'molecule': {
 'input_units_to_au': 1.8897261328856432,
 #'name': 'default',
}},
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
 #'units': 'Angstrom'
 }}


def test_efpefp1_b():
    """qchem-efp-sp"""

    pdict = molparse.init_psi4_molecule_from_any_string(_efpefp1_str)
    #print('\nRDICT')
    #pprint.pprint(_efpefp1_rdict)
    #print('\nWDICT')
    #pprint.pprint(pdict)
    assert(compare_dicts(_efpefp1_rdict, pdict, 6, sys._getframe().f_code.co_name + ': [1] str --> dict'))


@using_qcdb
def test_efpefp1_q():
    import qcdb
    pdict = molparse.init_psi4_molecule_from_any_string(_efpefp1_str)

    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
    assert(compare_values(0.0, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': [2] dict --> qcdbobj'))
    qdict = qmol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_efpefp1_rdict['molecule'])
    #print('\nQDICT')
    #pprint.pprint(qdict)
    assert(compare_dicts(_efpefp1_rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ': [3] qcdbobj --> dict'))


@using_pylibefp
def test_efpefp1_e():
    import pylibefp

    pdict = molparse.init_psi4_molecule_from_any_string(_efpefp1_str)

    emol = pylibefp.from_dict(pdict['libefp'])
    assert(compare_values(_efpefp1_nre, emol.nuclear_repulsion_energy(), 3, sys._getframe().f_code.co_name + ': [4] dict --> efpobj'))

    # too big
    #vdict = emol.to_viz_dict()
    #assert(compare_dicts(_efpefp1_vdict['libefp'], vdict, 4, sys._getframe().f_code.co_name + ': [5] efpobj --> viz dict'))

    edict = emol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_efpefp1_rdict['libefp'])
    #print('\nEDICT')
    #pprint.pprint(edict)
    assert(compare_dicts(_efpefp1_rdict['libefp'], edict, 4, sys._getframe().f_code.co_name + ': [6] efpobj --> dict'))


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
  'fragment_file': 'nh3'}],
 'molecule': {
 'input_units_to_au': 1.8897261328856432,
 #'name': 'default',
}},
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
 #'units': 'Angstrom'
 }}

_libefp1_vdict = {'libefp': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0, 0.0],
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
 #'units': 'Bohr'
 }}
# 'input_units_to_au': 1.8897261328856432,
# 'name': 'default',
# 'units': 'Angstrom'}}


def test_libefp1_b():

    pdict = molparse.init_psi4_molecule_from_any_string(_libefp1_str)
    #print('\nRDICT')
    #pprint.pprint(_libefp1_rdict)
    #print('\nWDICT')
    #pprint.pprint(pdict)
    assert(compare_dicts(_libefp1_rdict, pdict, 6, sys._getframe().f_code.co_name + ': [1] str --> dict'))


@using_qcdb
def test_libefp1_q():
    import qcdb
    pdict = molparse.init_psi4_molecule_from_any_string(_libefp1_str)

    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
    assert(compare_values(0.0, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': [2] dict --> qcdbobj'))

    qdict = qmol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_libefp1_rdict['molecule'])
    #print('\nQDICT')
    #pprint.pprint(qdict)
    assert(compare_dicts(_libefp1_rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ': [3] qcdbobj --> dict'))


@using_pylibefp
def test_libefp1_e():
    import pylibefp

    pdict = molparse.init_psi4_molecule_from_any_string(_libefp1_str)

    emol = pylibefp.from_dict(pdict['libefp'])
    assert(compare_values(_libefp1_nre, emol.nuclear_repulsion_energy(), 3, sys._getframe().f_code.co_name + ': [4] dict --> efpobj'))

    vdict = emol.to_viz_dict()
    #print('\nRDICT')
    #pprint.pprint(_libefp1_vdict['libefp'])
    #print('\nVDICT')
    #pprint.pprint(vdict)
    assert(compare_dicts(_libefp1_vdict['libefp'], vdict, 4, sys._getframe().f_code.co_name + ': [5] efpobj --> viz dict'))

    edict = emol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_libefp1_rdict['libefp'])
    #print('\nEDICT')
    #pprint.pprint(edict)
    assert(compare_dicts(filter_radian_range(_libefp1_rdict)['libefp'], edict, 4, sys._getframe().f_code.co_name + ': [6] efpobj --> dict'))


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
  'fragment_file': 'ammonia'}],
 'molecule': {'input_units_to_au': 1.8897261328856432}},
 'molecule': {
 'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0],
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
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': -1.423528806, 'y':  0.0,         'z': -0.897393146},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':  1.423528806, 'y':  0.0,         'z': -0.897393146}],
 'input_units_to_au': 1.8897261328856432,
 'name': 'default',
 'system_charge': 0.0,
 'system_multiplicity': 1,
 #'units': 'Angstrom'}}
 #'input_units_to_au': 1.0,
 #'units': 'Bohr'}}
 }}


_qmefp1_vdict = {'libefp': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0, 0.0, 0.0],
 'fragment_multiplicities': [1, 1, 1],
 'fragment_types': ['Real', 'Real', 'Real'],
 'fragments': [[0, 3], [3, 7], [7, 11]],
 'full_atoms': [{'charge': 8.0, 'ghosted': False, 'label': 'A01O1', 'mass': 15.99491461956, 'at_type': 'efpxyz', 'symbol': 'O', 'Z': 8.0,   'x': -4.043484801, 'y':  2.437066410, 'z': -1.822036143},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -5.043016189, 'y':  0.964404979, 'z': -2.163225699},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -2.519001456, 'y':  1.759578808, 'z': -1.114105301},
                {'charge': 7.0, 'ghosted': False, 'label': 'A01N1', 'mass': 14.00307400478, 'at_type': 'efpxyz', 'symbol': 'N', 'Z': 7.0,   'x':  1.866898241, 'y':  3.546656903, 'z':  5.389007602},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x':  3.189822098, 'y':  2.246050749, 'z':  5.849024393},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x':  2.756589676, 'y':  4.835894818, 'z':  4.293945633},
                {'charge': 1.0, 'ghosted': False, 'label': 'A04H4', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x':  1.431508879, 'y':  4.450724844, 'z':  7.015643157},
                {'charge': 7.0, 'ghosted': False, 'label': 'A01N1', 'mass': 14.00307400478, 'at_type': 'efpxyz', 'symbol': 'N', 'Z': 7.0,   'x': -7.800676093, 'y': -1.747354166, 'z': -2.426294971},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -8.868068505, 'y': -2.070327109, 'z': -3.978609780},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -6.787723797, 'y': -3.343360142, 'z': -2.144269196},
                {'charge': 1.0, 'ghosted': False, 'label': 'A04H4', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -9.023788093, 'y': -1.636376864, 'z': -0.961874594}],
 'input_units_to_au': 1.0},
# 'units': 'Bohr'},
'molecule': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0],
 'fragment_multiplicities': [1],
 'fragment_types': ['Real'],
 'fragments': [[0, 3]],
 'full_atoms': [{'charge': 8.0, 'ghosted': False, 'label': '0', 'mass': 15.99491461956, 'at_type': 'qmxyz', 'symbol': 'O', 'Z': 8.0,   'x':  0.0,           'y':  0.0,           'z':  0.224348286},
                {'charge': 1.0, 'ghosted': False, 'label': 'H_Aa', 'mass': 1.00782503207, 'at_type': 'qmxyz', 'symbol': 'H', 'Z': 1.0,   'x': -1.423528806, 'y':  0.0,           'z': -0.897393146},
                {'charge': 1.0, 'ghosted': False, 'label': 'H_Bb', 'mass': 1.00782503207, 'at_type': 'qmxyz', 'symbol': 'H', 'Z': 1.0,   'x':  1.423528806, 'y':  0.0,           'z': -0.897393146}],
 'input_units_to_au': 1.0}}
# 'units': 'Bohr'}}


def test_qmefp1_b():
    """qchem-qmefp-puream-sp"""

    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp1_str)
    #print('\nRDICT')
    #pprint.pprint(_qmefp1_rdict)
    #print('\nPDICT')
    #pprint.pprint(pdict)
    assert(compare_dicts(_qmefp1_rdict, pdict, 3, sys._getframe().f_code.co_name + ': [1] str --> dict'))


@using_qcdb
def test_qmefp1_q():
    import qcdb
    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp1_str)

    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
    qmol.print_out_in_bohr()
    assert(compare_values(_qmefp1_nre_qm, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': [2] dict --> qcdbobj'))

    qdict = qmol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_qmefp1_rdict['molecule'])
    #print('\nQDICT')
    #pprint.pprint(qdict)
    assert(compare_dicts(_qmefp1_rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ': [3] qcdbobj --> dict'))


@using_pylibefp
def test_qmefp1_e():
    import pylibefp

    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp1_str)

    emol = pylibefp.from_dict(pdict['libefp'])
    assert(compare_values(_qmefp1_nre_efp, emol.nuclear_repulsion_energy(), 4, sys._getframe().f_code.co_name + ': [4] dict --> efpobj'))
    # can't initialize qm charges in efp directly from molparse b/c mayn't be xyz
    emol.set_point_charges([at['Z'] for at in pdict['molecule']['full_atoms']],
                           [[at['x'], at['y'], at['z']] for at in pdict['molecule']['full_atoms']])

    assert(compare_values(_qmefp1_nre_qm, emol.nuclear_repulsion_energy(use_efp_frags=False, use_point_charges=True), 4, sys._getframe().f_code.co_name + ': [4b] qm geom'))
    assert(compare_values(_qmefp1_nre_tot, emol.nuclear_repulsion_energy(use_point_charges=True), 4, sys._getframe().f_code.co_name + ': [4c] qmefp geom'))

    vdict = emol.to_viz_dict()
    #print('\nRDICT')
    #pprint.pprint(_qmefp1_vdict['libefp'])
    #print('\nVDICT')
    #pprint.pprint(vdict)
    assert(compare_dicts(_qmefp1_vdict['libefp'], vdict, 4, sys._getframe().f_code.co_name + ': [5] efpobj --> viz dict'))

    # won't match b/c input is points, not xyzabc
    #edict = emol.to_dict()
    #assert(compare_dicts(filter_radian_range(_qmefp1_rdict)['libefp'], edict, 4, sys._getframe().f_code.co_name + ': [6] efpobj --> dict'))

_qmefp2_str = """
units au
0 1
O       0   0.0 0.
H100    2.0 0.0 0.0
H200    0.0 2.0 0.0
--
efp water -3.023563 8.881716,, 2.645618 -1.3 0.1 7
--
efp cl2    0.755891 -1.700754 -1.322809  2.3 1.6 -2.3
"""


_qmefp2_rdict = {'libefp': {'full_fragments':
 [{'coordinates_hint': [-3.023563, 8.881716, 2.645618, -1.3, 0.1, 7.0],
  'efp_type': 'xyzabc',
  'fragment_file': 'water'},
  {'coordinates_hint': [ 0.755891, -1.700754, -1.322809,  2.3, 1.6, -2.3],
  'efp_type': 'xyzabc',
  'fragment_file': 'cl2'}],
 'molecule': {'input_units_to_au': 1.0}},
 'molecule': {
 'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0],
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
                 'x':  0.0,         'y':  0.0,         'z':  0.0},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H100',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': 2.0, 'y':  0.0, 'z': 0.0},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H200',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'y': 2.0, 'x':  0.0, 'z': 0.0}],
 'input_units_to_au': 1.0,
 'name': 'default',
 'system_charge': 0.0,
 'system_multiplicity': 1,
 }}


_qmefp2_vdict = {'libefp': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0, 0.0],
 'fragment_multiplicities': [1, 1],
 'fragment_types': ['Real', 'Real'],
 'fragments': [[0, 3], [3, 5]],
 'full_atoms': [{'charge': 8.0, 'ghosted': False, 'label': 'A01O1', 'mass': 15.99491461956, 'at_type': 'efpxyz', 'symbol': 'O', 'Z': 8.0,   'x': -3.035639317257, 'y':  8.878363425377, 'z':  2.770530394361},
                {'charge': 1.0, 'ghosted': False, 'label': 'A02H2', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -2.149084890620, 'y': 10.094967315184, 'z':  1.761520373569},
                {'charge': 1.0, 'ghosted': False, 'label': 'A03H3', 'mass': 1.00782503207, 'at_type': 'efpxyz', 'symbol': 'H', 'Z': 1.0,   'x': -3.706381240250, 'y':  7.721672463301, 'z':  1.547265791286},
                {'charge': 17.0, 'ghosted': False, 'label': 'A01CL1', 'mass': 34.96885, 'at_type': 'efpxyz', 'symbol': 'Cl', 'Z': 17.0, 'x': -0.041833610179, 'y': -0.746973782608, 'z':  0.067493582349},
                {'charge': 17.0, 'ghosted': False, 'label': 'A02CL2', 'mass': 34.96885, 'at_type': 'efpxyz', 'symbol': 'Cl', 'Z': 17.0, 'x': 1.553615610179 , 'y': -2.654534217392, 'z': -2.713111582349}],
 'input_units_to_au': 1.0},
'molecule': {'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [0.0],
 'fragment_multiplicities': [1],
 'fragment_types': ['Real'],
 'fragments': [[0, 3]],
 'full_atoms': [{'charge': 8.0, 'ghosted': False, 'label': '0', 'mass': 15.99491461956, 'at_type': 'qmxyz', 'symbol': 'O', 'Z': 8.0,   'x':  0.0, 'y':  0.0, 'z':  0.0},
                {'charge': 1.0, 'ghosted': False, 'label': 'H100', 'mass': 1.00782503207, 'at_type': 'qmxyz', 'symbol': 'H', 'Z': 1.0,   'x': 2.0, 'y':  0.0, 'z': 0.0},
                {'charge': 1.0, 'ghosted': False, 'label': 'H200', 'mass': 1.00782503207, 'at_type': 'qmxyz', 'symbol': 'H', 'Z': 1.0,   'x': 0.0, 'y':  2.0, 'z': 0.0}],
 'input_units_to_au': 1.0}}


def test_qmefp2_b():

    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp2_str)
    #print('\nRDICT')
    #pprint.pprint(_qmefp2_rdict)
    #print('\nPDICT')
    #pprint.pprint(pdict)
    assert(compare_dicts(_qmefp2_rdict, pdict, 3, sys._getframe().f_code.co_name + ': [1] str --> dict'))


@using_qcdb
def test_qmefp2_q():
    import qcdb
    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp2_str)

    qmol = qcdb.Molecule.from_dict(pdict['molecule'])
    qmol.print_out_in_bohr()
    assert(compare_values(_qmefp2_nre_qm, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': [2] dict --> qcdbobj'))

    qdict = qmol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_qmefp2_rdict['molecule'])
    #print('\nQDICT')
    #pprint.pprint(qdict)
    assert(compare_dicts(_qmefp2_rdict['molecule'], qdict, 6, sys._getframe().f_code.co_name + ': [3] qcdbobj --> dict'))


@using_pylibefp
def test_qmefp2_e():
    import pylibefp

    pdict = molparse.init_psi4_molecule_from_any_string(_qmefp2_str)

    emol = pylibefp.from_dict(pdict['libefp'])
    assert(compare_values(_qmefp2_nre_efp, emol.nuclear_repulsion_energy(), 4, sys._getframe().f_code.co_name + ': [4] dict --> efpobj'))
    # can't initialize qm charges in efp directly from molparse b/c mayn't be xyz
    emol.set_point_charges([at['Z'] for at in pdict['molecule']['full_atoms']],
                           [[at['x'], at['y'], at['z']] for at in pdict['molecule']['full_atoms']])

    assert(compare_values(_qmefp2_nre_qm, emol.nuclear_repulsion_energy(use_efp_frags=False, use_point_charges=True), 4, sys._getframe().f_code.co_name + ': [4b] qm geom'))
    assert(compare_values(_qmefp2_nre_tot, emol.nuclear_repulsion_energy(use_point_charges=True), 4, sys._getframe().f_code.co_name + ': [4c] qmefp geom'))
    #print(emol.geometry_summary())

    vdict = emol.to_viz_dict()
    #print('\nRDICT')
    #pprint.pprint(_qmefp2_vdict['libefp'])
    #print('\nVDICT')
    #pprint.pprint(vdict)
    assert(compare_dicts(_qmefp2_vdict['libefp'], vdict, 4, sys._getframe().f_code.co_name + ': [5] efpobj --> viz dict'))

    edict = emol.to_dict()
    #print('\nRDICT')
    #pprint.pprint(_qmefp2_rdict['libefp'])
    #print('\nEDICT')
    #pprint.pprint(edict)
    assert(compare_dicts(filter_radian_range(_qmefp2_rdict)['libefp'], edict, 4, sys._getframe().f_code.co_name + ': [6] efpobj --> dict'))

