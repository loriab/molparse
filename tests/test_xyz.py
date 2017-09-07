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


_frag1_rdict = {'fragment_charges': [0, 0],
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
                 'x': 0.0,
                 'y': -0.667578,
                 'z': -2.124659},
                {'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x': 0.0,
                 'y': 0.667578,
                 'z': -2.124659},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': 0.923621,
                 'y': -1.232253,
                 'z': -2.126185},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': -0.923621,
                 'y': -1.232253,
                 'z': -2.126185},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': -0.923621,
                 'y': 1.232253,
                 'z': -2.126185},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': 0.923621,
                 'y': 1.232253,
                 'z': -2.126185},
                {'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x': 0.0,
                 'y': 0.0,
                 'z': 2.900503},
                {'Z': 6.0,
                 'charge': 6.0,
                 'ghosted': False,
                 'label': 'C',
                 'mass': 12,
                 'qm_type': 'qmcart',
                 'symbol': 'C',
                 'x': 0.0,
                 'y': 0.0,
                 'z': 1.69324},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': 0.0,
                 'y': 0.0,
                 'z': 0.627352},
                {'Z': 1.0,
                 'charge': 1.0,
                 'ghosted': False,
                 'label': 'H',
                 'mass': 1.00782503207,
                 'qm_type': 'qmcart',
                 'symbol': 'H',
                 'x': 0.0,
                 'y': 0.0,
                 'z': 3.963929}],
 'input_units_to_au': 1.8897261328856432,
 'name': 'dimer',
 'units': 'Angstrom'}


def test_frag1_b():
    """tu5"""

    wdict = molparse.init_psi4_molecule_from_any_string(_frag1_str, name='dimer')
    #print('RDICT', _frag1_rdict, '\n\n', 'WDICT', wdict, '\n')
    assert(compare_dicts(_frag1_rdict, wdict['molecule'], 6, sys._getframe().f_code.co_name + ': str --> dict'))


@using_qcdb
def test_frag1_q():
    import qcdb
    wdict = molparse.init_psi4_molecule_from_any_string(_frag1_str, name='dimer')

    emol = 85.189064196429

    qmol = qcdb.Molecule.from_dict(wdict['molecule'])
    assert(compare_values(emol, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': dict --> qcdbobj'))

    qdict = qmol.to_dict()
    #print('RDICT', _frag1_rdict, '\n\n', 'QDICT', qdict, '\n')
    assert(compare_dicts(_frag1_rdict, qdict, 6, sys._getframe().f_code.co_name + ': qcdbobj --> dict'))


_frag2_str = """

He 0 0 0
--
@He 0 0 4
--
He 0 4 0
units au
"""


_frag2_rdict = {'fragment_charges': [0, 0, 0],
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
 'name': 'default',
 'units': 'Bohr'}


def test_frag2_b():
    """nbody-he-cluster"""

    wdict = molparse.init_psi4_molecule_from_any_string(_frag2_str)
    #print('RDICT', _frag2_rdict, '\n\n', 'WDICT', wdict, '\n')
    assert(compare_dicts(_frag2_rdict, wdict['molecule'], 6, sys._getframe().f_code.co_name + ': str --> dict'))


@using_qcdb
def test_frag2_q():
    import qcdb
    wdict = molparse.init_psi4_molecule_from_any_string(_frag2_str)

    emol = 1.000

    qmol = qcdb.Molecule.from_dict(wdict['molecule'])
    assert(compare_values(emol, qmol.nuclear_repulsion_energy(), 6, sys._getframe().f_code.co_name + ': dict --> qcdbobj'))

    qdict = qmol.to_dict()
    #print('RDICT', _frag2_rdict, '\n\n', 'QDICT', qdict, '\n')
    assert(compare_dicts(_frag2_rdict, qdict, 6, sys._getframe().f_code.co_name + ': qcdbobj --> dict'))


## <<<  https://stackoverflow.com/a/18860653
#def dict_compare(d1, d2):
#    d1_keys = set(d1.keys())
#    d2_keys = set(d2.keys())
#    intersect_keys = d1_keys.intersection(d2_keys)
#    added = d1_keys - d2_keys
#    removed = d2_keys - d1_keys
#    modified = {o : (d1[o], d2[o]) for o in intersect_keys if d1[o] != d2[o]}
#    same = set(o for o in intersect_keys if d1[o] == d2[o])
#    return added, removed, modified, same
## https://stackoverflow.com/a/18860653 >>>
#
#
#def dsame(ref, test):
#    """Compares dictionary *test* to *ref* by keys and values.
#
#    Uses dict_compare.
#    Returns True if no keys added or removed in *test* with respect to
#    *ref* and if all the values are identical, else False.
#
#    """
#    added, removed, modified, same = dict_compare(test, ref)
#    if added or removed or modified:
#        if added:
#            print('  Added:', *added)
#        if removed:
#            print('  Removed:', *removed)
#        if modified:
#            print('  Modified:')
#            for item in modified:
#                try:
#                    if isinstance(test[item][0], dict):
#                        for idx, ll in enumerate(test[item]):
#                            print('  {}: '.format(idx), end='')
#                            if dsame(ref[item][idx], test[item][idx]):
#                                print('same')
#                except TypeError:
#                #else:
#                    print('    ', item, ':', ref[item], '-->', test[item])
#        return False
#    else:
#        return True


if __name__ == "__main__":
    #test_frag1()
    #test_frag2()

    test_frag1_b()
    test_frag1_q()
    test_frag2_b()
    test_frag2_q()

