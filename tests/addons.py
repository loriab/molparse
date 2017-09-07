import pytest

def _plugin_import(plug):
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


using_qcdb = pytest.mark.skipif(_plugin_import('qcdb') is False,
                                reason='molparse not detecting package qcdb. Install package if necessary and and to envvar PYTHONPATH')
using_pylibefp = pytest.mark.skipif(_plugin_import('pylibefp') is False,
                                reason='molparse not detecting package pylibefp. Install package if necessary and and to envvar PYTHONPATH')
using_psi4 = pytest.mark.skipif(_plugin_import('psi4') is False,
                                reason='molparse not detecting package psi4. Install package if necessary and and to envvar PYTHONPATH')


#@using_cfour
#def test_cfour():
#    """cfour/sp-rhf-ccsd_t_"""

