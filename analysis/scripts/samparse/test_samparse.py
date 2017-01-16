# Unit Test Module for Testing samparse!
# Nathan Lubock

import json
from parse_lib import *
from samparse import updateUncertain

#===============================================================================

def readJSON(path, tup='answer'):
    """
    A quick wrapper funtion to read in JSON file specifing the parameters for
    testing. Since the ouput of of samparse's functions are lists of tuples
    and JSON does not support them, we must edit them. We will edit them
    in-place for now.
    Assumes:
        JSON file is a list of dicts
        Contains one key:value pair that is a list of tups
    """
    with open(path, 'r') as f:
        store = json.load(f)
    for x in store:
        # Protect against None case when updateRead will skip certain reads
        if x[tup] is not None:
            x[tup] = [tuple(y) for y in x[tup]]
    return store

#===============================================================================

def pytest_generate_tests(metafunc):
    # called once per each test function
    for funcargs in metafunc.cls.params[metafunc.function.__name__]:
        # schedule a new test function run with applied **funcargs
        metafunc.addcall(funcargs=funcargs)

#-------------------------------------------------------------------------------

class TestRetrieval:
    """
    Set of tests for any function related to finding different classes of
    errors. Note that the numbering here should be with resepect to the read
    itself, not the reference sequence. Since each parameter contains the
    settings for all tests in general, we must filter out the relevant parts
    of the answer.
    """
    store = readJSON('./unitTest/retrieval_Params.json')
    params = {
            'test_getInserts':[{key:x[key] for key in ('name', 'diff', 'cigar', 'answer')} for x in store],
            'test_getMM':[{key:x[key] for key in ('name', 'diff', 'cigar', 'md', 'answer')} for x in store],
            'test_getDels':[{key:x[key] for key in ('name', 'diff', 'cigar', 'md', 'answer')} for x in store]
            }

    # NOTE: answer refers to actual output of samparse, we want output of funs!
    def test_getInserts(self, name, diff, cigar, answer):
        out = [x for x in answer if x[1] in ('I', 'S')]
        assert out == list(getInserts(diff, cigar))

    def test_getMM(self, name, diff, cigar, md, answer):
        out = [x for x in answer if x[1] == 'M']
        assert out == list(getMM(diff, cigar, md))

    def test_getDels(self, name, diff, cigar, md, answer):
        out = [x for x in answer if x[1] in ('D', 'P')]
        my_diff = add_s(diff, cigar)
        assert out == list(getDels(cigar, md))

#-------------------------------------------------------------------------------

# class TestUpdatePos:
#     """
#     Set of tests relevant to the helper function updatePos
#     """
#     store = readJSON('./unitTest/updatePos_Params.json')
#     params = {'test_updatePos':store}
#
#     def test_updatePos(self, name, raw_list, diff, ref, left, answer):
#         assert answer == updatePos(raw_list, diff, ref, left)

#-------------------------------------------------------------------------------

class TestUpdateRead:
    """
    Set of tests relevant to the final function updateRead. See:
    https://holgerkrekel.net/2009/05/13/parametrizing-python-tests-generalized/
    for explination on how this works.
    """
    updateRead_Params = readJSON('./unitTest/updateRead_Params.json')
    updateUncertain_Params = readJSON('./unitTest/updateUncertain_Params.json')

    # we can simply read in store for now...
    params = {'test_updateRead':updateRead_Params,
              'test_updateUncertain':updateUncertain_Params}

    def test_updateRead(self, name, ref, diff, cigar, md, left, answer):
        # protect against None case
        # if answer is not None:
        #     out = (name, answer)
        # else:
        #     out = answer
        assert answer == updateRead(name, ref, diff, cigar, md)

    def test_updateUncertain(self, name, canon, ref, answer):
        # no need to test None case, since we clean them in humpty
        # convert canon to a list of tuples since json's can't store them
        # (answer is automatically converted by readJSON)
        canon = [tuple(x) for x in canon]
        assert answer == updateUncertain(canon, ref)
