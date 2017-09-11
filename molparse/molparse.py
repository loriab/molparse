#"""Module with utility functions that act on molecule objects."""
from __future__ import absolute_import
import re
#import math
#
#from psi4 import core
#from psi4.driver.p4util import constants
#from psi4.driver.inputparser import process_pubchem_command, pubchemre
#from psi4.driver import qcdb
#from psi4.driver.qcdb import periodictable
from . import constants
from . import periodictable

NUMBER = r'((?:[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?)|(?:[-+]?\d+\.\d*(?:[DdEe][-+]?\d+)?)|(?:[-+]?\d+(?:[DdEe][-+]?\d+)?))'
SEP = r'[\t ,]+'


    #molecule = init_psi4_molecule_from_any_string(geom, name=name)


def init_psi4_molecule_from_any_string(mol_str, **kwargs):



#    #if not has_efp:

    # Figure out how and if we will parse the Molecule adata
    mname = kwargs.pop("name", "default")
    dtype = kwargs.pop("dtype", "psi4").lower()
    if mol_str is not None:

        print('<<< QMMOL', mol_str, '>>>')
        # << 1 >> str -- discard comments
        mol_str = filter_comments(mol_str)

        if dtype == 'psi4':

            # << 2 >> dict[m] -- seed Psi4 minimal and default fields
            mol_init = {'name': mname,
                        #'units': 'Angstrom',
                        #'input_units_to_au': 1.0 / constants.bohr2angstroms,
                        #'units': 'Bohr',
                        'input_units_to_au': 1.0 / constants.bohr2angstroms,
                        'fragments': [],
                        'fragment_types': [],
                        'fragment_charges': [],
                        'fragment_multiplicities': []}

            # << 3 >>  str-->dict[m] -- process units, com, orient, symm
            mol_str, univ = filter_universals(mol_str)
            mol_init.update(univ)

            # << 4 >>  str-->dict[e] -- process efp
            mol_str, efp_init = filter_libefp(mol_str, confer=mol_init)
            if efp_init:
                print('<<< core.EFP INTO', efp_init, '>>>')
                # GOOD! core.efp_init()
                # GOOD! efp = core.get_active_efp()
                # GOOD! efp.construct_from_pydict(efp_init)
                # << 5 >>  dict[e]-->dict[m] --- tie qm & efp axes
                mol_init['fix_com'] = True
                mol_init['fix_orientation'] = True
                mol_init['fix_symmetry'] = 'c1'

            mol_str, mints_init = filter_mints(mol_str, confer=mol_init)
            mol_init.update(mints_init)

            mol_init = reconcile_cgmp(mol_init)
            print('\nINTO core.Molecule.from_dict <<<', mol_init, '>>>\n')
            return {'molecule': mol_init, 'libefp': efp_init}
#            molecule = core.Molecule.from_dict(mol_init)
#            molecule.update_geometry()
#            molecule.print_out()
#            return molecule
#
            #has_efp, geom = filter_libefp(geom)

        else:
            raise KeyError("Molecule: dtype of %s not recognized.")
    else:
        # In case a user wants to build one themselves
        pass


def reconcile_cgmp(pydict):
    def apply_default(lst, default):
        return [default if (c is None) else c for c in lst]

    def reconcile_charges(chg, fchg):
        """expects None as placeholder for not specified"""

        naive_fchg = apply_default(fchg, 0.0)

        if chg is None:
            # Sys: None, Frag: [-1, 1] --> Sys: 0, Frag: [-1, 1]
            # Sys: None, Frag: [None, -2] --> Sys: -2, Frag: [0, -2]
            return sum(naive_fchg), naive_fchg
        else:
            if None in fchg:
                # Sys: 0, Frag: [None, 1] --> Sys: 0, Frag: [-1, 1]
                # Sys: 2, Frag: [1, None, None] --> Sys: 2, Frag: [1, 1, 0]
                missing_frag_chg = chg - sum(naive_fchg)
                first = fchg.index(None)
                tmp = fchg[:]
                tmp[first] = missing_frag_chg
                return chg, apply_default(tmp, 0.0)
            else:
                if chg == sum(fchg):
                    # Sys: 0, Frag: [-1, 1] --> Sys: 0, Frag: [-1, 1]
                    return chg, fchg
                else:
                    # Sys: 0, Frag: [-2, 1] --> Irreconcilable!
                    raise ValidationError('System charge: {} not reconcilable with fragment charges: {}'.format(chg, fchg))

    def reconcile_multiplicities(mult, fmult):

        def high_spin_sum(mult_list):
            mm = 1
            for m in mult_list:
                mm += m - 1
            return mm

        naive_fmult = apply_default(fmult, 1)

        if mult is None:
            return high_spin_sum(naive_fmult), naive_fmult
        else:
            if None in fmult:
                missing_frag_mult = mult - high_spin_sum(naive_fmult) + 1
                first = fmult.index(None)
                tmp = fmult[:]
                tmp[first] = missing_frag_mult
                return mult, apply_default(tmp, 1)
            else:
                if mult == high_spin_sum(fmult):
                    return mult, fmult
                else:
                    # TODO: could be ok, could be not
                    raise ValidationError('Spin Arithmetic')

    chg = pydict.get('system_charge', None)
    mult = pydict.get('system_multiplicity', None)
    fchg = pydict['fragment_charges']
    fmult = pydict['fragment_multiplicities']
    # TODO look at the atoms and see if chgmult appropriate

    if (chg is not None) or (len(fchg) == 1 and fchg[0] is not None):
        chg_specified = True
    else:
        chg_specified = False

    if (mult is not None) or (len(fmult) == 1 and fmult[0] is not None):
        mult_specified = True
    else:
        mult_specified = False

    chg, fchg = reconcile_charges(pydict.get('system_charge'), pydict['fragment_charges'])
    mult, fmult = reconcile_multiplicities(pydict.get('system_multiplicity'), pydict['fragment_multiplicities'])

    if chg_specified:
        pydict['system_charge'] = chg
    pydict['fragment_charges'] = fchg
    if mult_specified:
        pydict['system_multiplicity'] = mult
    pydict['fragment_multiplicities'] = fmult

    return pydict


def filter_pubchem(mol_str):
    pass
    pubchemerror = re.compile(r'^\s*PubchemError\s*$', re.IGNORECASE)
    pubcheminput = re.compile(r'^\s*PubchemInput\s*$', re.IGNORECASE)


def filter_comments(mol_str):
    comment = re.compile(r'(^|[^\\])#.*')
    mol_str = re.sub(comment, '', mol_str)
    return mol_str


def filter_universals(mol_str):  # needed?, seed=None):
    """Process multiline string *mol_str* for fix_ and unit markers,
    returning a tuple of unprocessed *mol_str* and a dictionary of
    processed fields.

    fix_com
    fix_orientation
    fix_symmetry
    input_units_to_au
    units

    """
    com = re.compile(r'\A(no_com|nocom)\Z', re.IGNORECASE)
    orient = re.compile(r'\A(no_reorient|noreorient)\Z', re.IGNORECASE)
    bohrang = re.compile(r'\Aunits?[\s=]+((?P<ubohr>(bohr|au|a.u.))|(?P<uang>(ang|angstrom)))\Z', re.IGNORECASE)
    symmetry = re.compile(r'\Asymmetry[\s=]+(?P<pg>\w+)\Z', re.IGNORECASE)

    def process_com(matchobj):
        univ['fix_com'] = True
        return ''

    def process_orient(matchobj):
        univ['fix_orientation'] = True
        return ''

    def process_bohrang(matchobj):
        if matchobj.group('uang'):
            #univ['units'] = 'Angstrom'
            univ['input_units_to_au'] = 1.0 / constants.bohr2angstroms
        elif matchobj.group('ubohr'):
            #univ['units'] = 'Bohr'
            univ['input_units_to_au'] = 1.0
        return ''

    def process_symmetry(matchobj):
        univ['fix_symmetry'] = matchobj.group('pg').lower()
        return ''

    reconstitute = []
    univ = {}
    # below needed?
    #if seed:
    #    univ = seed

    for line in mol_str.split('\n'):
        line = re.sub(com, process_com, line.strip())
        line = re.sub(orient, process_orient, line)
        line = re.sub(bohrang, process_bohrang, line)
        line = re.sub(symmetry, process_symmetry, line)
        if line:
            reconstitute.append(line)

    return '\n'.join(reconstitute), univ


#def filter_mints(mol_str, seed=None):
def filter_mints(mol_str, confer=None):
    """confer is read-only"""

    fragment_re = re.compile(r'^\s*--\s*$', re.MULTILINE)
    ATOM = r'(?:(?P<gh1>@)|(?P<gh2>Gh\())?(?P<label>(?P<symbol>[A-Z]{1,3})(?:(_\w+)|(\d+))?)(?(gh2)\))(?:@(?P<mass>\d+\.\d+))?'
    CHG = r'(?P<chg>-?\d+)'
    MULT = r'(?P<mult>\d+)'
    cgmp = re.compile(r'\A' + CHG + SEP + MULT + r'\Z')
    atom_cart = re.compile(r'\A' + ATOM + SEP + NUMBER + SEP + NUMBER + SEP + NUMBER + r'\Z', re.IGNORECASE)
    atom_zmat1 = re.compile(r'\A' + ATOM + r'\Z', re.IGNORECASE)
    atom_zmat2 = re.compile(r'\A' + ATOM + SEP + r'(\d+)' + SEP + NUMBER + r'\Z', re.IGNORECASE)
    atom_zmat3 = re.compile(r'\A' + ATOM + SEP + r'(\d+)' + SEP + NUMBER
                                         + SEP + r'(\d+)' + SEP + NUMBER + r'\Z', re.IGNORECASE)
    atom_zmat4 = re.compile(r'\A' + ATOM + SEP + r'(\d+)' + SEP + NUMBER
                                         + SEP + r'(\d+)' + SEP + NUMBER
                                         + SEP + r'(\d+)' + SEP + NUMBER + r'\Z', re.IGNORECASE)
    #frag = re.compile(r'^\s*--\s*$')
    #variable = re.compile(r'^\s*(\w+)\s*=\s*(-?\d+\.\d+|-?\d+\.|-?\.\d+|-?\d+|tda)\s*$', re.IGNORECASE)
    #ghost = re.compile(r'@(.*)|Gh\((.*)\)', re.IGNORECASE)

    def process_system_cgmp(matchobj):
        mints_init['system_charge'] = int(matchobj.group('chg'))
        mints_init['system_multiplicity'] = int(matchobj.group('mult'))
        return ''

    def filter_fragment(mol_str):

        def process_fragment_cgmp(matchobj):
            mints_init['fragment_charges'].append(float(matchobj.group('chg')))
            mints_init['fragment_multiplicities'].append(int(matchobj.group('mult')))
            return ''

        def process_atom_spec(matchobj):
            atom_symbol = matchobj.group('symbol').upper()
            if atom_symbol not in periodictable.el2z:
                raise ValidationError("""Illegal atom symbol in geometry specification: {}""".format(atom_symbol))

            if matchobj.group('gh1') or matchobj.group('gh2'):
                atom_Z = 0.0
                atom_charge = 0.0
                atom_ghosted = True
            else:
                atom_Z = float(periodictable.el2z[atom_symbol])
                atom_charge = float(atom_Z)
                atom_ghosted = False

            if matchobj.group('mass'):
                atom_mass = float(matchobj.group('mass'))
            else:
                atom_mass = periodictable.el2mass[atom_symbol]

            atom_init = {}
            atom_init['Z'] = atom_Z
            atom_init['charge'] = atom_charge
            atom_init['ghosted'] = atom_ghosted
            atom_init['symbol'] = atom_symbol.capitalize()
            atom_init['label'] = matchobj.group('label').capitalize()
            atom_init['mass'] = atom_mass
            return atom_init

        def process_atom_cart(matchobj):

            atom_init = process_atom_spec(matchobj)
            atom_init['qm_type'] = 'qmcart'
            atom_init['x'] = float(matchobj.group(8)) * input_units_to_au
            atom_init['y'] = float(matchobj.group(9)) * input_units_to_au
            atom_init['z'] = float(matchobj.group(10)) * input_units_to_au

            mints_init['full_atoms'].append(atom_init)
            return ''

        #moldict['fragment_charges'] = self.fragment_charges
        #moldict['fragment_multiplicities'] = self.fragment_multiplicities

        #frag_reconstitute = []
        #start_atom = len(mints_init["full_atoms"])

        #for iln, line in enumerate(mol_str.split('\n')):
        #    line = line.strip()
        #    if iln == 0:
        #        line = re.sub(cgmp, process_fragment_cgmp, line)
        #    line = re.sub(atom_cart, process_atom_cart, line)
        #    #line = re.sub(atom_zmat1, process_atom_zmat1, line)
        #    #line = re.sub(atom_zmat2, process_atom_zmat2, line)
        #    #line = re.sub(atom_zmat3, process_atom_zmat3, line)
        #    #line = re.sub(atom_zmat4, process_atom_zmat4, line)
        #    if line:
        #        frag_reconstitute.append(line)

        #end_atom = len(mints_init["full_atoms"])
        #if end_atom > 0:
        #    #mints_init['fragments'].append([start_atom, end_atom + 1])
        #    mints_init['fragments'].append([start_atom, end_atom])
        #    mints_init['fragment_types'].append('Real')
        #    if len(mints_init['fragment_charges']) < len(mints_init['fragments']):
        #        mints_init['fragment_charges'].append(0)
        #        mints_init['fragment_multiplicities'].append(1)

        frag_reconstitute = []
        start_atom = len(mints_init["full_atoms"])

        for iln, line in enumerate(mol_str.split('\n')):
            line = re.sub(atom_cart, process_atom_cart, line.strip())
            #line = re.sub(atom_zmat1, process_atom_zmat1, line)
            #line = re.sub(atom_zmat2, process_atom_zmat2, line)
            #line = re.sub(atom_zmat3, process_atom_zmat3, line)
            #line = re.sub(atom_zmat4, process_atom_zmat4, line)
            if line:
                frag_reconstitute.append(line)

        end_atom = len(mints_init["full_atoms"])
        if (end_atom - start_atom) > 0:
            if re.sub(cgmp, process_fragment_cgmp, mol_str.split('\n')[0].strip()) != '':
                mints_init['fragment_charges'].append(None)
                mints_init['fragment_multiplicities'].append(None)
            #mints_init['fragments'].append([start_atom, end_atom + 1])
            mints_init['fragments'].append([start_atom, end_atom])
            mints_init['fragment_types'].append('Real')

        return '\n'.join(frag_reconstitute), mints_init

    reconstitute = []
    mints_init = {}
    mints_init['fragments'] = []
    mints_init['fragment_types'] = []
    mints_init['fragment_charges'] = []
    mints_init['fragment_multiplicities'] = []
    mints_init['full_atoms'] = []
    #if seed:
    #    mints_init = seed

    if confer and ('input_units_to_au' in confer):
        input_units_to_au = confer['input_units_to_au']
    else:
        input_units_to_au = 1.0

    # handle `--`-demarcated blocks
    for ifr, frag in enumerate(re.split(fragment_re, mol_str)):
        frag = frag.strip()
        if ifr == 0:
            frag = re.sub(cgmp, process_system_cgmp, frag)
        frag, mints_init = filter_fragment(frag)
        if frag:
            reconstitute.append(frag)

    return '\n--\n'.join(reconstitute), mints_init


def filter_libefp(mol_str, confer=None):
    """confer is read-only"""

    ENDL = r'[\t ,]*$'

    fragment = re.compile(r'^\s*--\s*$', re.MULTILINE)
    efpxyzabc = re.compile(
        r'\A' + r'efp' + SEP + r'(\w+)' +
        SEP + NUMBER + SEP + NUMBER + SEP + NUMBER +
        SEP + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL + r'\Z',
        re.IGNORECASE)
    efppoints = re.compile(
        r'\A' + r'efp' + r'\s+' + r'(\w+)' + ENDL +
        r'[\s,]*' + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL +
        r'[\s,]*' + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL +
        r'[\s,]*' + NUMBER + SEP + NUMBER + SEP + NUMBER + ENDL + r'\Z',
        re.IGNORECASE | re.MULTILINE)

    def process_efpxyzabc(matchobj):
        efp_frags.append({'efp_type': 'xyzabc',
        # NOTE: imposing case on file
                         'fragment_file': matchobj.group(1).lower(),
                         'coordinates_hint': [float(matchobj.group(2)) * input_units_to_au,
                                              float(matchobj.group(3)) * input_units_to_au,
                                              float(matchobj.group(4)) * input_units_to_au,
                                              float(matchobj.group(5)),
                                              float(matchobj.group(6)),
                                              float(matchobj.group(7))]})
        return ''

    def process_efppoints(matchobj):
        efp_frags.append({'efp_type': 'points',
                         'fragment_file': matchobj.group(1).lower(),
                         'coordinates_hint': [float(matchobj.group(2)) * input_units_to_au,
                                              float(matchobj.group(3)) * input_units_to_au,
                                              float(matchobj.group(4)) * input_units_to_au,
                                              float(matchobj.group(5)) * input_units_to_au,
                                              float(matchobj.group(6)) * input_units_to_au,
                                              float(matchobj.group(7)) * input_units_to_au,
                                              float(matchobj.group(8)) * input_units_to_au,
                                              float(matchobj.group(9)) * input_units_to_au,
                                              float(matchobj.group(10)) * input_units_to_au]})
        return ''

    reconstitute = []
    efp_init = {}
    efp_frags = []  # list of dicts, one for each efp fragment in mol_str

    # << 4.1 >>  findsldkjfsl
    # NOTE: applying libefp default of AU
    if confer and ('input_units_to_au' in confer):
        input_units_to_au = confer['input_units_to_au']
    else:
        input_units_to_au = 1.0

    # handle `--`-demarcated blocks
    for frag in re.split(fragment, mol_str):
        frag = re.sub(efpxyzabc, process_efpxyzabc, frag.strip())
        frag = re.sub(efppoints, process_efppoints, frag)
        if frag:
            reconstitute.append(frag)
    if efp_frags:
        efp_init['full_fragments'] = efp_frags
        efp_init['molecule'] = {'input_units_to_au': input_units_to_au}

    return '\n--\n'.join(reconstitute), efp_init


