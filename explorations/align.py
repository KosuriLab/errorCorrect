#!/data/home/nlubock/miniconda2/bin/python

"""
align.py - a wrapper around UTA Align for testing alignment accuracy
Nathan Lubock
"""
import itertools
import argparse
import multiprocessing
import sys
import csv
from uta_align.align.algorithms import align, needleman_wunsch_altshul_erikson
from signal import signal, SIGPIPE, SIG_DFL
from functools import partial

# catch broken pipe errors to allow ex) python parse.py foo bar | head
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


def fastaReader(path):
    """
    Read in a fasta file lazily and return a generator of the name and sequence

    Parameters:
    -----------
    path :: str
        path to fasta file

    Returns:
    --------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        seq :: str
            Sequence taken from the fasta file

    Dependencies:
    -------------
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
    with open(path, 'r') as fp:
        # ditch the boolean (x[0]) and just keep the header/seq grouping
        fa_iter = (x[1] for x in itertools.groupby(fp, lambda line: line[0] == ">"))
        for header in fa_iter:
            # drop the ">"
            name = header.next()[1:].strip()
            # join all sequence lines to one by iterating until the next group.
            seq = "".join(s.strip() for s in fa_iter.next())
            yield name, seq


#-------------------------------------------------------------------------------


def fastqReader(path):
    """
    Read in a fastq file laizy and return a generator of the name and sequence

    Parameters:
    -----------
    path :: str
        path to fasta file

    Returns:
    --------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        seq :: str
            Sequence taken from the fasta file

    Dependencies:
    -------------
    grouper - see itertools recepies
    """
    with open(path, 'r') as fp:
        group_gen = grouper(fp, 4)
        for record in group_gen:
            name = record[0].split(' ')[0][1:].strip()
            seq = record[1].strip()
            yield name, seq


#===============================================================================

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

    return (name, align_out.gapped_alignment()[0], align_out.gapped_alignment()[1])


#-------------------------------------------------------------------------------


if __name__ == '__main__':
    # create an argument parser a nicer CLI
    parser = argparse.ArgumentParser(
            description='Process input and output as a csv)')
    parser.add_argument('ref_fasta',
            help='path to a *.fasta file of the reference sequence')
    parser.add_argument('query',
            help='file of queries')
    parser.add_argument('-a', '--algo',
            dest='algo', nargs='?', default='global', metavar='foo',
            choices=['global', 'local', 'glocal', 'local_global', 'altshul'],
            help='alignment algorithm to use (default=global, choices=global, local, glocal, local_global, altshul). Please refer to source code for more info on each')
    parser.add_argument('-m', '--match', dest='match', type=int,
            default=10, metavar='N',
            help='match score (default=10)')
    parser.add_argument('-x', '--mismatch', dest='mismatch', type=int,
            default=-9, metavar='N',
            help='mismatch penalty (default=-9)')
    parser.add_argument('-o', '--gap-open', dest='open_pen', type=int,
            default=-15, metavar='N',
            help='gap open penalty (default=-15)')
    parser.add_argument('-e', '--gap-extend', dest='extend_pen', type=int,
            default=-6, metavar='N',
            help='gap extend penalty (default=-6)')
    parser.add_argument('-p', '--proc', dest='proc', type=int, default=1,
            metavar='N', choices=range(1, multiprocessing.cpu_count()),
            help='number of processors (default=1, max={})'.format(multiprocessing.cpu_count()))
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
    ref = list(next(fastaReader(args.ref_fasta)))[1]

    # check for fasta/fastq query (hackish and brittle...)
    if args.query.split('.')[-1] == 'fasta':
        query_gen = fastaReader(args.query)
    else:
        query_gen = fastqReader(args.query)

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
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['Name', 'Ref', 'Seq'])
    writer.writerows(errs)
