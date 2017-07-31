"""
parse.py - a simple CLI for getting all errors in an alingment.

TODO:
-----
    * Make sure input fasta/fastq's are valid
    * Parse file extensions to automatically read fasta/fastq
    * Exit gracefully on KeyboardInterrupt
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import re
import itertools
import argparse
import multiprocessing
import sys
import csv
from uta_align.align.algorithms import align, needleman_wunsch_altshul_erikson
from signal import signal, SIGPIPE, SIG_DFL

# catch broken pipe errors to allow eg) python parse.py foo bar | head
# see: http://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================
# recepies/helper functions


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)

def flatten(listOfLists):
    "Flatten one level of nesting"
    return itertools.chain.from_iterable(listOfLists)


#===============================================================================
# I/O Functions

def fasta_reader(fasta):
    """
    Read in a fasta file lazily and return a generator of the name and sequence

    Parameters:
    -----------
    fasta :: FileType
        opened file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        read :: str
            Sequence taken from the fasta file

    Requires:
    ---------
    itertools

    Example:
    --------
    itertools.groupby takes a key function and groups all items into a list
    until that key changes. We can key on lines beginning with >, then grab
    every line until the next record in the fasta. This makes our method robust
    to some fasta formats that have forced line breaks at given characters.

    foo = '>ABC>DEF>GHI'
    [(k, list(g)) for k,g in itertools.groupby(foo, lambda x: x == '>')]
    --> [(True, ['>']), (False, ['A', 'B', 'C']), (True, ['>']), ... ]

    Note:
    -----
    Adapted from: https://www.biostars.org/p/710/#1412
    """
    # ditch the boolean (x[0]) and just keep the header/seq grouping
    fa_iter = (x[1] for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read


def fastq_reader(fastq):
    """
    Read in a fastq file laizy and return a generator of the name and sequence

    Parameters:
    -----------
    fastq :: FileType
        opened fastq file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        seq :: str
            Sequence taken from the fasta file

    Dependencies:
    -------------
    grouper - see itertools recepies
    """
    group_gen = grouper(fastq, 4)
    for record in group_gen:
        # drop the @ before the name and any text after a whitespace
        name = record[0].split(' ')[0][1:].strip()
        seq = record[1].strip()
        yield name, seq




#===============================================================================
# Parsing Functions


def updateIns(err_list):
    """
    Cumulatively renumber the positions of the errors caused by insertions
    padding the reference sequence.

    Parameters:
    -----------
    err_list :: [(Pos, Type, Diff), ...]
        Sorted list of tuples corresponding to our errors

    Returns:
    --------
    mod_list :: [(Pos, Type, Diff), ...]

    Example:
    --------
    Update the length of all differences by compensating for insertions. We can
    think of the algorithm as such:

        0, M, ...  = 0  - 0
        2, I, A    = 2  - 0
        7, M, ...  = 7  - A
        9, I, B    = 9  - A
        12, M, ... = 12 - A - B
        ...

    Basically, we subtract the length of an insertion from the positions of
    every error downstream of it. If there are more than one insertion, we
    subtract the cumulative length of the two insertions.

    """
    cum_len = 0
    mod_list = []
    for i in err_list:
        mod_list.append((i[0] - cum_len, i[1], i[2]))
        if i[1] == 'I' or i[1] == 'S':
            cum_len = cum_len + len(i[2])
    return mod_list


#-------------------------------------------------------------------------------


def getErrs(name, ref, query):
    """
    Get the relevant errors from two input strings

    Parameters
    ----------
    name :: str
        Name of read. Unused, only here for the output data structure
    ref :: str
        Reference sequence of the alignment. Must use -'s for padding
        insertions. Must have actual sequence in [A-Za-z]
    query :: str
        Query of the alignment. Must use -'s for padding deletions.
        Matching characters to reference must be replaced with non [A-Z].
        Mismatches must be in [A-Z]

    Returns
    -------
    [(Name, Pos, Type, Diff), ...]
        Name :: Str
            Name of read
        Pos :: Int
            Position of Error (0 indexed)
        Type :: Str
            M - mismatch
            D - deletion
            P - multiPle deletions
            I - insertion
            S - multiple inSertions
        Diff :: Str
            M - original base and what it is mismatched to
            D|P - base(s) that were deleted
            I|S - base(s) that are inserted

    Dependencies:
    -------------
    itertools
    re
    """
    # use re to find indel positions, and MM will be what's left behind
    # remember ins pad the ref sequence, while dels pad the query
    ins_pos = [m.span() for m in re.finditer('-+', ref)]
    del_pos = [m.span() for m in re.finditer('-+', query)]
    letter_pos = [m.span() for m in re.finditer('[A-Z]', query)]

    mm_pos = set(flatten(itertools.starmap(xrange, letter_pos))) - \
             set(flatten(itertools.starmap(xrange, ins_pos)))

    # grab the actual errors from the ref/query and handle multiple ins/dels
    ins_gen = ((x[0], query[slice(*x)]) for x in ins_pos)
    del_gen = ((x[0], ref[slice(*x)]) for x in del_pos)

    ins_list = [(x[0], 'I', x[1]) if len(x[1]) == 1 else (x[0], 'S', x[1]) for x in ins_gen]
    del_list = [(x[0], 'D', x[1]) if len(x[1]) == 1 else (x[0], 'P', x[1]) for x in del_gen]
    mm_list = [(x, 'M', ref[x] + query[x]) for x in mm_pos]

    # correct error position caused by insertions then convert to 1-based index
    combined = sorted(itertools.chain(mm_list, del_list, ins_list), key = lambda x: x[0])
    if ins_list == []:
        return [(name, x[0] + 1, x[1], x[2]) for x in combined]
    else:
        return [(name, x[0] + 1, x[1], x[2]) for x in updateIns(combined)]


#-------------------------------------------------------------------------------


def updateUncertain(update_read, my_ref):
    """
    Take single base insertions and deletions in repeat regions and distribute
    their "count" over the regio. For example, if we deleted an 'A' from the
    sequence 'TAAAG', we would not know which 'A' that was. We will find these
    cases and note the error at all valid positions.

    Paramerters:
    ------------
    update_read :: [(Name, Pos, Type, Diff), ...]
        List of every error for a given read packaged into a list of tuples
    my_ref :: str
        Reference sequence to find repeated regions

    Returns:
    --------
    [(Name, Pos, Type, Diff), ...]
        List of tuples sorted on Pos. Note that for uncertain types, Diff will
        be converted to $N$X where:
            $N - Length of uncertain region
            $X - Original Error
        and Type will be converted to U$T where $T was the original type

    Example:
    --------
    ref:  TAAAG
    read: T-AAG
    update_read = [('foo', 2, 'D', 'A')]
    return = [('foo', 2, 'UD', '3A'),
              ('foo', 3, 'UD', '3A'),
              ('foo', 4, 'UD', '3A')]

    Dependencies:
    -------------
    re
    itertools
    flatten
    """
    # get positions of repeats in the ref_seq
    # (\w) catures a 'word' char, \1 back-reference, {1,} reps \1 >= 1 times
    match = re.finditer(r'(\w)\1{1,}', my_ref)
    rep_pos = [xrange(*m.span()) for m in match]

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
    rng_pack = [zip(rng, itertools.repeat(len(rng), len(rng))) for rng in rng_gen]

    # create new tuples with uncertain tag in first pass
    # add length of unncertain on the second pass
    # return numbering to 1-based
    new_tups = ((tup[0], x[0] + 1, 'U' + tup[2], str(x[1]) + tup[3])
            for i, tup in enumerate(uncertain) for x in rng_pack[i])

    # combine the new tuples with the good ones and sort on the index
    return sorted(itertools.chain(good, new_tups), key=lambda x: x[1])


#-------------------------------------------------------------------------------


def parseWrap(pack):
    """
    Transform alignment output into a form amenable for analysis. Note, we have
    to package everything together to get around multiprocessing's limitations.
    Specifically, pool.map requires the function to be pickle-able.
    Unfortunately, we cannot pickle lambdas or partially applied functions, so
    we must pass all parameters as one giant tuple.

    Parameters:
    -----------
    pack :: (ref, algo, {match, mismatch, open, extend}, (name, query))
        ref :: str
            Reference sequence for alignment
        algo :: str
            What aligner to use (global, local, glocal, local_global, altshul).
            See Notes for more details
        match :: int
            Match score (default 10)
        mismatch :: int
            Mismatch pentalty (default -9)
        open :: int
            Gap open penalty (default -15)
        extend :: int
            Gap extend penalty (default -6)
        name :: str
            Name of read
        query :: srt
            Sequence to align

    Returns:
    --------
    [(Name, Pos, Type, Diff), ...]
        List of tuples pertaining to each error in the query relative to the
        reference sequence. Will be ordered by position.

    Notes:
    ------
    (Taken from uta_align source)
    Global - Needleman-Wunsch Gotoh based global alignment
    Local - Smith-Waterman Gotoh based local alignment
    Glocal - "similar to global alignment with no penalty for leading or
        trailing gaps, provided that the maximum scoring alignments that spans
        from the start of the reference or the query sequence and read the end
        of the reference or the query sequence"
    Local_global - "similar to global alignment with no penalty for leading gaps
        and to glocal alignments except the alignment must span to the end of
        both the reference and query alignments"
    Altshul - a variant of the standard Needleman-Wunsch global alignment. See:
        Altschul SF, Erickson BW.  "Optimal sequence alignment using affine gap
        costs." Bull Math Biol.  1986;48(5-6):603-16.
    """
    ref, algo, score_dict, read_pack = pack
    name, query = read_pack

    if algo == 'altshul':
        align_out = needleman_wunsch_altshul_erikson(ref, query, **score_dict)
    else:
        align_out = align(ref, query, algo, **score_dict)

    errs = getErrs(name, align_out.gapped_alignment()[0], align_out.gapped_alignment()[1])
    uncert = updateUncertain(errs, ref)
    return uncert


#-------------------------------------------------------------------------------


if __name__ == '__main__':
    # create an argument parser a nicer CLI
    parser = argparse.ArgumentParser(
            description='Process input and output as a csv')
    parser.add_argument('ref_fasta',
            type=argparse.FileType('r'),
            help='path to a *.fastA file of the reference sequence')
    parser.add_argument('query',
            type=argparse.FileType('r'),
            default=sys.stdin,
            nargs='?',
            help='path to a *.fastQ file of the reads (or stdin if none)')
    parser.add_argument('-p',
            '--proc',
            dest='proc',
            type=int,
            default=1,
            metavar='N',
            choices=range(1, multiprocessing.cpu_count() + 1),
            help='number of processors (default=1, max={})'.format(multiprocessing.cpu_count()))
    parser.add_argument('-a',
            '--algo',
            dest='algo',
            nargs='?',
            default='global',
            metavar='foo',
            choices=['global', 'local', 'glocal', 'local_global', 'altshul'],
            help='alignment algorithm to use (default=global, choices=global,' +
            'local, glocal, local_global, altshul). Please refer to source' +
            'code for more info on each')
    parser.add_argument('-m',
            '--match',
            dest='match',
            type=int,
            default=10,
            metavar='N',
            help='match score (default=10)')
    parser.add_argument('-x',
            '--mismatch',
            dest='mismatch',
            type=int,
            default=-9,
            metavar='N',
            help='mismatch penalty (default=-9)')
    parser.add_argument('-o',
            '--gap-open',
            dest='open_pen',
            type=int,
            default=-15,
            metavar='N',
            help='gap open penalty (default=-15)')
    parser.add_argument('-e',
            '--gap-extend',
            dest='extend_pen',
            type=int,
            default=-6,
            metavar='N',
            help='gap extend penalty (default=-6)')
    parser.add_argument('-fa',
            '--fasta',
            dest='fasta',
            action='store_true',
            help='process query as fasta file')
    args = parser.parse_args()

    # Make sure that the score values are all within a valid range
    if args.match <= args.mismatch:
        raise ValueError('Match score must be greater than mismatch penalty!')
    if args.match <= args.open_pen:
        raise ValueError('Match score must be greater than gap open penalty!')
    if args.match <= args.extend_pen:
        raise ValueError('Match score must be greater than gap extend penalty!')
    if args.open_pen > args.extend_pen:
        raise ValueError('Gap open penalty must be less than or equal to gap extend penalty!')

    # process input
    ref = next(fasta_reader(args.ref_fasta))[1]

    # check for fasta/fastq query (hackish and brittle...)
    if args.fasta:
        query_gen = fasta_reader(args.query)
    else:
        query_gen = fastq_reader(args.query)

    # pack everything into a neat generator
    score_dict = {'match_score':args.match, 'mismatch_score':args.mismatch,
            'gap_open_score':args.open_pen, 'gap_extend_score':args.extend_pen}
    pack = ((ref, args.algo, score_dict, x) for x in query_gen)

    # set up pool and distribute
    pool = multiprocessing.Pool(args.proc)
    errs = pool.map(parseWrap, pack, chunksize=10000)
    pool.close()
    pool.join()

    # output csv to stdout
    tidy_df = flatten(x for x in errs if x != [])
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['Name', 'Pos', 'Type', 'Diff'])
    writer.writerows(tidy_df)
