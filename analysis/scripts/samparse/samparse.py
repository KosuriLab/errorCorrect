#!/data/home/nlubock/miniconda2/envs/py35/bin/python
"""
samparse.py - CLI to run the parser
Nathan Lubock
"""

import pandas as pd
import itertools as it
import argparse
import re
import sys
import csv
from functools import partial
from parse_lib import updateRead
from signal import signal, SIGPIPE, SIG_DFL

# catch broken pipe errors to allow eg) python samparse.py foo bar | head
# see: http://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)


#===============================================================================
# Aliases/Helper Functions


def flatten(iterable):
    """
    Shorten itertools obnoxiously long itertools.chain.from_iterable
    """
    return it.chain.from_iterable(iterable)

#-------------------------------------------------------------------------------

def fastaReader(path):
    """
    Get the reference sequence from the fasta file specified by the user. This
    solution is somewhat cheeky and rigid. However, since we are dealing with a
    known fasta file, with a single entry, we only need the second line!
    """
    with open(path, 'r') as f:
        next(f)
        return f.readline().rstrip()

#-------------------------------------------------------------------------------

def samReader(path, aligner, rows=None):
    """
    Read a .sam file into a PANDAS data frame into memory. Note columns are
    hard-coded! Manually assigns types... may not be needed.
    Make lazy in the future...
    """
    md_dict = {'bowtie':17, 'bbmap':13, 'needle':13}
    md_pos = md_dict[aligner]

    raw_read = pd.read_table(path, usecols=[0, 3, 4, 5, 9, 10, md_pos],
            names=['Name', 'Left', 'MapQ', 'CIGAR', 'Read', 'BaseQ', 'MD'],
            dtype={'Name': str, 'Left': int, 'MapQ': str, 'CIGAR': str,
                'Read': str, 'BaseQ': str, 'MD': str},
            skiprows=range(0, 3), header=None, nrows=rows)

    # filter out md garbage in md col
    raw_read['MD'] = raw_read.ix[0:, 'MD'].str[5:]
    return raw_read

#===============================================================================
# I/O functions

def updateUncertain(update_read, my_ref):
    """
    Take single base insertions and deletions in repeat regions and distribute
    their "count" over the region. When there are deletions or insertions in
    repeat regions, Bowtie will default to the left-most position. However, it
    is equally valid that the error occurs at any one of the positions in the
    repeat. For example, if we deleted an 'A' from the sequence 'TAAAG', we
    would not know which 'A' that was. We will find these cases and note the
    error at all valid positions.

    Input:
            read_update - list of form [(Name, Pos, Type, Diff), ...]
    Output:
            list of tups with (Name, Pos, Type, Diff)
    """
    # get positions of repeats in the ref_seq
    # (\w) catures a 'word' char, \1 back-reference, {1,} reps \1 >= 1 times
    match = re.finditer(r'(\w)\1{1,}', my_ref)
    rep_pos = [range(m.start(), m.end()) for m in match]

    # make sure the indexing matches that of the reference sequence
    # (updateRead will output 1 based, but we need 0 based) we will do this
    # in-place for now
    update_read = [(x[0], x[1] - 1, x[2], x[3]) for x in update_read]

    # split update read into tuples that are in repeat regions or are not
    # make sure good numbering returns to a 1-index
    rep_set = set(flatten(rep_pos))
    uncertain = [x for x in update_read
            if x[1] in rep_set and x[2] in set(['D', 'I'])]
    good = ((x[0], x[1] + 1, x[2], x[3]) for x in update_read
            if x[1] not in rep_set or x[2] not in set(['D', 'I']))

    # find the repeat region for each tuple in uncertain and package each
    # position and the length of the range so we can do everything in one pass
    # remember list comp loops go outer then inner
    rng_gen = (rng for rng in rep_pos for tup in uncertain if tup[1] in rng)
    rng_pack = [zip(rng, it.repeat(len(rng), len(rng))) for rng in rng_gen]

    # create new tuples with uncertain tag in first pass
    # add length of unncertain on the second pass
    # return numbering to 1-based
    new_tups = ((tup[0], x[0] + 1, 'U' + tup[2], str(x[1]) + tup[3])
            for i, tup in enumerate(uncertain) for x in rng_pack[i])

    # combine the new tuples with the good ones and sort on the index
    return list(sorted(it.chain(good, new_tups), key=lambda x: x[1]))


#===============================================================================
# Main body


if __name__ == '__main__':
    # create an argument parser a nicer CLI
    parser = argparse.ArgumentParser(
            description='Process *.sam file for easier analysis. Will output a tidy csv to stdout')
    parser.add_argument('sams',
            help='path to the *.sam file')
    parser.add_argument('ref_fasta',
            help='path to a *.fasta file of the reference sequence')
    parser.add_argument('aligner', nargs='?', default='bowtie',
            choices=['bowtie', 'bbmap', 'needle'],
            help='Aligner used in pipeline (hard coded relevant columns in sam file for now. Defaults to bowtie')
    parser.add_argument('-u', '--uncertain', dest='canon', action='store_false',
            help='distribute errors over repeat regions?')
    args = parser.parse_args()

    if args.canon == False:
        print('NOTE: Will distribute errors over regions with repeated bases',
                file=sys.stderr)

    # attempt to open files
    ref = fastaReader(args.ref_fasta)
    read_df = samReader(args.sams, args.aligner)

    # itertuple returns (_, Name, Left, MapQ, Cigar, Read, BaseQ, MD)
    # updateRead(my_name, my_ref, my_diff, my_cigar, my_md):
    pack = ((f[1], args.ref_fasta, f[5], f[4], f[7]) for f in read_df.itertuples())

    raw_reads = (updateRead(*x) for x in pack)
    clean_reads = (x for x in raw_reads if x is not None)

    if args.canon == False:
        tidy_df = (updateUncertain(x, ref) for x in clean_reads)
    else:
        tidy_df = clean_reads

    # output to stdout
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['Name', 'Pos', 'Type', 'Diff'])
    writer.writerows(flatten(tidy_df))
