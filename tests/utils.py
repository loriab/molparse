import sys
import math
import pprint
import deepdiff


def _success(label):
    """Function to print a '*label*...PASSED' line to screen.
    Used by :py:func:`util.compare_values` family when functions pass.

    """
    msg = '\t{0:.<66}PASSED'.format(label)
    print(msg)
    sys.stdout.flush()


# Test functions
def compare_values(expected, computed, digits, label, exitonfail=True):
    """Function to compare two values. Prints :py:func:`util.success`
    when value *computed* matches value *expected* to number of *digits*
    (or to *digits* itself when *digits* < 1 e.g. digits=0.04). Performs
    a system exit on failure unless *exitonfail* False, in which case
    returns error message. Used in input files in the test suite.

    """
    if digits > 1:
        thresh = 10 ** -digits
        message = ("\t%s: computed value (%.*f) does not match (%.*f) to %d digits." % (label, digits+1, computed, digits+1, expected, digits))
    else:
        thresh = digits
        message = ("\t%s: computed value (%f) does not match (%f) to %f digits." % (label, computed, expected, digits))
    if abs(expected - computed) > thresh:
        print(message)
        if exitonfail:
            return False
    if math.isnan(computed):
        print(message)
        print("\tprobably because the computed value is nan.")
        if exitonfail:
            return False
    _success(label)
    return True


def compare_integers(expected, computed, label, exitonfail=True):
    """Function to compare two integers. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if (expected != computed):
        print("\t%s: computed value (%d) does not match (%d)." % (label, computed, expected))
        return False
    _success(label)
    return True


def compare_dicts(expected, computed, tol, label):
    """Compares dictionaries *computed* to *expected* using DeepDiff

    Float comparisons made to *tol* significant decimal places.

    Note that a clean DeepDiff returns {}, which evaluates to False, hence the compare_integers.

    """
    ans = deepdiff.DeepDiff(expected, computed, significant_digits=tol, verbose_level=2)
    clean = not bool(ans)
    if not clean:
        pprint.pprint(ans)
    return compare_integers(True, clean, label)


## <<<  https://stackoverflow.com/a/18860653
#def dict_compare(d1, d2):
#    d1_keys = set(d1.keys())
#    d2_keys = set(d2.keys())
#    intersect_keys = d1_keys.intersection(d2_keys)
#    added = d1_keys - d2_keys
#    removed = d2_keys - d1_keys
#    #modified = {o : (d1[o], d2[o]) for o in intersect_keys if d1[o] != d2[o]}
#    #same = set(o for o in intersect_keys if d1[o] == d2[o])
#    modified = {}
#    same = []
#    for o in intersect_keys:
#        if isinstance(d1[o], float):
#            if abs(d1[o] - d2[o]) < 1e-5:
#                same.append(o)
#            else:
#                modified[o] = (d1[o], d2[o])
#        else:
#            if d1[o] == d2[o]:
#                same.append(o)
#            else:
#                modified[o] = (d1[o], d2[o])
#    same = set(same)
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
#            truly_modified = []
#            for item in modified:
#                #print('  ', item)
#                try:
#                    if isinstance(test[item], dict):
#                        if dsame(ref[item], test[item]):
#                            print('same2')
#                        else:
#                            return False
#                    elif isinstance(test[item][0], dict):
#                        cleans = []
#                        for idx, ll in enumerate(test[item]):
#                            print('  {}: '.format(idx), end='')
#                            if dsame(ref[item][idx], test[item][idx]):
#                                print('same')
#                                cleans.append(True)
#                        if len(cleans) != len(test[item]):
#                            truly_modified.append(item)
#                    else:
#                        print('    ', item, ':', ref[item], '-->', test[item])
#                        truly_modified.append(item)
#                except TypeError:
#                    print('    ', item, ':', ref[item], '-->', test[item])
#                    truly_modified.append(item)
#        if truly_modified:
#            return False
#        else:
#            return True
#    else:
#        return True

