"""
parse_lib.py - a suite of functions to analyze sequences specified in a sam file
Nathan Lubock
"""

import operator
import re
import itertools as it


#=================================================================================
# Aliases/Helper Functions


def scanl(f, base, l):
    """
    Pythonic implementation of Haskell's scanl. Useful for cumsums.
    """
    yield base
    for x in l:
        base = f(base, x)
        yield base

#-------------------------------------------------------------------------------

def flatten(iterable):
    """
    Shorten itertools obnoxiously long itertools.chain.from_iterable
    """
    return it.chain.from_iterable(iterable)

#-------------------------------------------------------------------------------

def seqDiffs(my_seq):
    """
    Get all differences in a read with some regex.
    Idea:
        regex [A-Z]|_+ = all A-Z or greedy match _
        ex) ===A==T___==G -> [(3,'A'), (5,'T'), (6,'___'), (12,G)]
    """
    return ((m.start(0),  m.group()) for m in re.finditer('[A-Z]|_+', my_seq))


#===============================================================================
# Cigar string parsing functions


def splitTag(my_tag):
    """
    Split the CIGAR string into various operations with some regex. Assumes the
    tag is alphanumeric!
    Idea:
        regex (\d+) = repeats of digits, (\D) = [^0-9] (not numbers)
        (\D+) allows for matching of deletion characters
        ex) 30M1I4M -> [('30','M'), ('1','I'), ('4','M')]
    Note:
        WILL LEAVE OFF STRAGGLING NUMERICS
        ex) 30M1I4M4 -> [('30','M'), ('1','I'), ('4','M')]
    """
    my_split = re.findall(r'(\d+)(\D+)', my_tag)
    return ((int(x[0]), x[1]) for x in my_split)

#-------------------------------------------------------------------------------

def cumUpdateCigar(my_cigar):
    """
    Cumulatively update the numbering of a cigar string starting at 0
    Input:
        my_cigar - literal CIGAR string
    Output:
        List of tuples in form of splitTag
    """
    # fancy tuple unpacking with zip
    # zip(*$VAR) is the inverse of zip
    my_pos, my_str = zip(*splitTag(my_cigar))

    # cumsum the positions of the cigar strings
    # note this will return a list of n+1 since it starts at 0
    new_pos = scanl(operator.add, 0, my_pos)

    # note zip ignores the added length of new_pos
    return zip(new_pos, my_str)

#-------------------------------------------------------------------------------

def getDiffLocs(my_cigar, kind):
    """
    Get the locations of your difference of choice (m,d,i).
    Note the range return is left inclusive [0,89) = [0,88] like python ranges
    Input:
        my_cigar - literal CIGAR string
        kind - string of what type you would like to get the location of
    Output:
        list of tuples with start and end of difference
    """
    # the cigar string already provides the length of differences
    # cumulatively add the lengths of the diffs to get their start positions
    # coerce cum_start to list since you will use it more than once
    diff_len = (x[0] for x in splitTag(my_cigar) if x[1] == kind)
    cum_start = [x[0] for x in cumUpdateCigar(my_cigar) if x[1] == kind]

    # recall map(f(x,y), (x1, x2), (y1, y2)) = (f(x1, y1), f(x2, y2))
    cum_end = map(operator.add, diff_len, cum_start)
    return zip(cum_start, cum_end)


#=============================================================================
# Retrieval functions
# return the position of a given class of errors


def getInserts(my_diff, my_cigar):
    """
    Get the location and makeup of all insertions in your read. Simply finds the
    locations of the inserts from the CIGAR string, then subsets the padded diff
    string to get the actual errors.

    Note - requires my_diff to be padded with _'s. Also remember that
    getDiffLocs returns tuples of beginning and end.
    """
    my_ranges = list(getDiffLocs(my_cigar, 'I'))
    my_ins = (my_diff[x:y] for x,y in my_ranges)
    return ((x[0], 'I', y) if len(y) == 1 else (x[0], 'S', y)
            for x,y in zip(my_ranges, my_ins))

#-------------------------------------------------------------------------------

def getMM(my_diff, my_cigar, my_md):
    """
    Intersect the positions of errors in the CIGAR String with the Diff String
    to find mismatches since Sam1.3 does not differentiate between matches and
    mismatches. Format the output as:

        [(pos, 'M', '$Orig$Miss),...]

    Use index from seqDiffs to be consistent by keeping everything in terms of
    the reference sequence.
    """
    # use starmap to get ranges since getDiffLocs only returns (start, end)
    # flatten the ranges into a set for quick lookup (list would work fine too)
    my_range = set(flatten(it.starmap(range, getDiffLocs(my_cigar, 'M'))))
    diff_pack = seqDiffs(my_diff)
    my_filter = [x for x in diff_pack if x[0] in my_range]

    # this protects against the bowtie error case where my_filter will be empty
    if my_filter == []:
        return my_filter
    else:
        idx, bases = zip(*my_filter)

    # filter out anything thats not a MM from the md tag
    orig_bases = (x[1] for x in splitTag(my_md) if x[1][0] != '^')
    mm = (x + y for x,y in zip(orig_bases, bases))
    return ((x, 'M', y) for x,y in zip(idx, mm))

#-------------------------------------------------------------------------------

def getMM_sam14(my_diff, my_cigar, my_md):
    """
    Compatibility function for Sam1.4 files. Since mismatches and matches are
    differentiated by 'X' vs '=', we can use the CIGAR string to get the
    locations, then use the diff string/md tag to get the mismatcha and original
    base. Output in form:

        [(pos, 'M', '$Orig$Miss),...]
    """
    mm_range = (x for x in getDiffLocs(my_cigar, 'X'))
    orig_bases = (x[1] for x in splitTag(my_md) if x[1][0] != '^')
    return [(x[0], 'M', y + my_diff[slice(*x)]) for x,y in zip(mm_range, orig_bases)]

#-------------------------------------------------------------------------------

def getDels(my_cigar, my_md):
    """
    Find the positions of the deletions with the CIGAR string, then use the MD
    tag to get the bases that were deleted
    """
    # we only need the first position of the getDiffLocs range
    # x[1][1:] cleans '^' from dels
    del_loc = (x[0] for x in getDiffLocs(my_cigar, 'D'))
    del_type = (x[1][1:] for x in splitTag(my_md) if x[1][0] == '^')
    return ((x, 'D', y) if len(y) == 1 else (x, 'P', y)
            for x,y in zip(del_loc, del_type))

#=============================================================================
# Update functions
# Used to get numbering of read seq to match the ref seq


def add_s(my_diff, my_cigar):
    """
    Insert a '_' to the read sequence for every del to ensure the length of the
    read stays 100
    """
    # since getDiffLocs outputs ((start, end),...), star map range over output
    del_pos = it.chain(*it.starmap(range, getDiffLocs(my_cigar, 'D')))
    split_diff = list(my_diff)
    for f in del_pos:
        split_diff.insert(f, '_')
    return ''.join(split_diff)

#-------------------------------------------------------------------------------

def updateIns(read_update, my_diff):
    """
    Update the length of all differences by compensating for insertions. We can
    think of the algorithm as such:

        0, L, ...  = 0
        2, I, A    = 2
        7, M, ...  = 7  - A
        9, I, B    = 9  - A
        12, M, ... = 12 - (A + B)
        ...

    Basically, we subtract the length of an insertion from the positions of
    every error downstream of it. If there are more than one insertion, we
    subtract the cumulative length of the two insertions.

    NOTE - read_update must be a list since it is called twice
    """
    # unpack the original read_update into vecs for downstream processing
    read_pos, kind, chars = zip(*read_update)

    # filter the position and length of each insertion into their own vectors
    # add one to the position for downstream ranges
    ins_pos, ins_len = zip(*((x[0] + 1, len(x[2])) for x in read_update
        if x[1] in ['I', 'S']))

    # use accumulate instead of scanl since we dont care about subtracting zero
    cum_len = list(it.accumulate(ins_len))

    # append the length of our read so we can zip ranges together
    # note the trailing comma ensures that this is a tuple concatenation
    # forgive the mutable state
    ins_pos = ins_pos + (len(my_diff) + 1,)
    ins_range = [range(x,y) for x,y in zip(ins_pos, ins_pos[1:])]

    out_pos = list(read_pos)
    for i in range(0, len(read_pos)):
        for j in zip(ins_range, cum_len):
            if out_pos[i] in j[0]:
                out_pos[i] = out_pos[i] - j[1]

    return zip(out_pos, kind, chars)


#===============================================================================
# Helper functions


def isMM(my_md):
    """
    Quick mismatch checker that looks for characters in the MD Tag.
    Further checks that the first base is a nt since splitTag greedy matches
    all non-digits (ex 50^ATG -> (50, ^ATG))
    """
    if True in [c in my_md for c in ['A','T','C','G']]:
        _, base = zip(*splitTag(my_md))
        return (x[0] in ['A', 'T', 'G', 'C'] for x in base)
    else:
        return [False]

#===============================================================================

def updateRead(my_name, my_ref, my_diff, my_cigar, my_md):
    """
    Find all of the errors and output them into an easy to parse key:value
    format.
    """
    # Flag for updating positions based on insert length downstream
    was_ins = False

    # Bail out if there is no match at all
    if my_cigar == '*':
        return None
    if 'D' in my_cigar:
        # major side effects will insert _'s into current diff_string
        my_diff = add_s(my_diff, my_cigar)
        del_update = getDels(my_cigar, my_md)
    else:
        del_update = []
    if 'I' in my_cigar:
        was_ins = True
        ins_update = getInserts(my_diff, my_cigar)
    else:
        ins_update = []
    if True in isMM(my_md):
        # Sam1.4 compatibility check
        if 'X' in my_cigar:
            mm_update = getMM_sam14(my_diff, my_cigar, my_md)
        else:
            mm_update = getMM(my_diff, my_cigar, my_md)
    else:
        mm_update = []

    # combine everything into a big list
    read_updates = list(it.chain(ins_update, mm_update, del_update))

    # exit early for perfects
    if read_updates == []:
        return None
    else:
        sorted_updates = sorted(read_updates, key = lambda x: x[0])
        if was_ins:
            final_updates = list(updateIns(sorted_updates, my_diff))
        else:
            final_updates = sorted_updates
        return [(my_name, x[0] + 1, x[1], x[2]) for x in final_updates]
