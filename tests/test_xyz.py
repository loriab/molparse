import molparse
from utils import *
from addons import *


_frag1_str = """
0 1
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
0 1
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
units angstrom
"""


_frag1_rdict = {'fragment_charges': [0.0, 0.0],
 'fragment_multiplicities': [1, 1],
 'fragment_types': ['Real', 'Real'],
 'fragments': [[0, 6], [6, 10]],
 'full_atoms': [{'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x':   0.000000000000, 'y':    -1.261539592340,
                 #'z':    -4.021632563049},  # com
                 'z':  -4.015023635771},  # no_com
                 #'x': 0.0,
                 #'y': -0.667578,
                 #'z': -2.124659},
                {'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x':  -0.000000000000, 'y':     1.261539592340,
                 #'z':    -4.021632563049},  # com
                 'z':  -4.015023635771},  # no_com
                 #'x': 0.0,
                 #'y': 0.667578,
                 #'z': -2.124659},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':   1.745390740582, 'y':    -2.328620696427,
                 #'z':    -4.024516285128},  # com
                 'z':  -4.017907357849},  # no_com
                 #'x': 0.923621,
                 #'y': -1.232253,
                 #'z': -2.126185},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':  -1.745390740582, 'y':    -2.328620696427,
                 #'z':    -4.024516285128},  # com
                 'z':  -4.017907357849},  # no_com
                 #'x': -0.923621,
                 #'y': -1.232253,
                 #'z': -2.126185},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':  -1.745390740582, 'y':     2.328620696427,
                 #'z':    -4.024516285128},  # com
                 'z':  -4.017907357849},  # no_com
                 #'x': -0.923621,
                 #'y': 1.232253,
                 #'z': -2.126185},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':   1.745390740582, 'y':     2.328620696427,
                 #'z':    -4.024516285128},  # com
                 'z':  -4.017907357849},  # no_com
                 #'x': 0.923621,
                 #'y': 1.232253,
                 #'z': -2.126185},
                {'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x':  -0.000000000000, 'y':     0.000000000000,
                 #'z':     5.474547390335},  # com
                 'z':   5.481156317613},  # no_com
                 #'x': 0.0,
                 #'y': 0.0,
                 #'z': 2.900503},
                {'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x':  -0.000000000000, 'y':     0.000000000000,
                 #'z':     3.193150949969},  # com
                 'z':   3.199759877247},  # no_com
                 #'x': 0.0,
                 #'y': 0.0,
                 #'z': 1.69324},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':  -0.000000000000, 'y':     0.000000000000,
                 #'z':     1.178914541639},  # com
                 'z':   1.185523468918},  # no_com
                 #'x': 0.0,
                 #'y': 0.0,
                 #'z': 0.627352},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x':  -0.000000000000, 'y':     0.000000000000,
                 #'z':     7.484131292925}],  # com
                 'z':   7.490740220203}],  # no_com
                 #'x': 0.0,
                 #'y': 0.0,
                 #'z': 3.963929}],
 'input_units_to_au': 1.8897261328856432,
 'name': 'dimer'}
# 'units': 'Angstrom'}


def test_frag1_b():
    """tu5"""

    wdict = molparse.init_psi4_molecule_from_any_string(_frag1_str, name='dimer')
    #print('\nRDICT')
    #pprint.pprint(_frag1_rdict)
    #print('\nWDICT')
    #pprint.pprint(wdict['molecule'])
    assert(compare_dicts(_frag1_rdict, wdict['molecule'], 6, sys._getframe().f_code.co_name + ': str --> dict'))


@using_qcdb
def test_frag1_q():
    import qcdb
    wdict = molparse.init_psi4_molecule_from_any_string(_frag1_str, name='dimer')

    emol = 85.189064196429

    qmol = qcdb.Molecule.from_dict(wdict['molecule'])
    assert(compare_values(emol, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': dict --> qcdbobj'))

    qdict = qmol.to_dict()
    assert(compare_dicts(_frag1_rdict, qdict, 6, sys._getframe().f_code.co_name + ': qcdbobj --> dict'))


_frag2_str = """

He 0 0 0
--
@He 0 0 4
--
He 0 4 0
units au
"""


_frag2_rdict = {'fragment_charges': [0.0, 0.0, 0.0],
 'fragment_multiplicities': [1, 1, 1],
 'fragment_types': ['Real', 'Real', 'Real'],
 'fragments': [[0, 1], [1, 2], [2, 3]],
 'full_atoms': [{'Z': 2.0,
                 'charge': 2.0,
                 'ghosted': False,
                 'label': 'He',
                 'mass': 4.00260325415,
                 'qm_type': 'qmcart',
                 'symbol': 'He',
                 'x': 0.0,
                 'y': 0.0,
                 'z': 0.0},
                {'Z': 0.0,
                 'charge': 0.0,
                 'ghosted': True,
                 'label': 'He',
                 'mass': 4.00260325415,
                 'qm_type': 'qmcart',
                 'symbol': 'He',
                 'x': 0.0,
                 'y': 0.0,
                 'z': 4.0},
                {'Z': 2.0,
                 'charge': 2.0,
                 'ghosted': False,
                 'label': 'He',
                 'mass': 4.00260325415,
                 'qm_type': 'qmcart',
                 'symbol': 'He',
                 'x': 0.0,
                 'y': 4.0,
                 'z': 0.0}],
 'input_units_to_au': 1.0,
 'name': 'default'}
# 'units': 'Bohr'}


def test_frag2_b():
    """nbody-he-cluster"""

    wdict = molparse.init_psi4_molecule_from_any_string(_frag2_str)
    assert(compare_dicts(_frag2_rdict, wdict['molecule'], 6, sys._getframe().f_code.co_name + ': str --> dict'))


@using_qcdb
def test_frag2_q():
    import qcdb
    wdict = molparse.init_psi4_molecule_from_any_string(_frag2_str)

    emol = 1.000

    qmol = qcdb.Molecule.from_dict(wdict['molecule'])
    assert(compare_values(emol, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': dict --> qcdbobj'))

    qdict = qmol.to_dict()
    assert(compare_dicts(_frag2_rdict, qdict, 6, sys._getframe().f_code.co_name + ': qcdbobj --> dict'))

