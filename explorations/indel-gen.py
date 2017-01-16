#!/data/home/nlubock/miniconda2/envs/py35/bin/python
"""
Insert or delete arbitrarily long N-mers at a given position
Nathan Lubock
June, 2016

TODO - Convert to command line utility
"""
import itertools
import csv
import sys
import argparse
from signal import signal, SIGPIPE, SIG_DFL

# catch broken pipe errors to allow ex) python parse.py foo bar | head
# see: http://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================

def flatten(x):
    return(itertools.chain.from_iterable(x))

def insert_nmers_at(pos, nmer, seq):
    """
    Insert nmers into your DNA sequence of interest. Note that the
    position of the insert is right indexed. E.G:
        ref: A_TCG
        seq: AGTCG (Insert G at 1)
        ref: _ATGC
        seq: GATGC (Insert G at 0)
    Input:
        seq  - string
        nmer - positive integer
        pos  - int in [0, len(seq)]
    Output:
        List of strings containing all nmers inserted at pos
    Reqs:
        itertools (imported as it)
        flatten alias for itertools.chain.from_iterable
    """
    nmer_list = list(itertools.product('ATGC', repeat=nmer))

    # combine all nmers at pos together
    out_seqs = (''.join(flatten(
        (itertools.islice(seq, 0, pos), x, itertools.islice(seq, pos, None))
        )) for x in nmer_list)

    # name format Ins_(nmer)_(@pos)
    out_names = ("Ins_{}_{}".format(''.join(x), pos) for x in nmer_list)

    # generate sam tag
    # out_list will automatically place ins at last possible base
    if pos == 0:
        out_sam = "{}I{}M".format(nmer, len(seq))
    elif pos >= len(seq):
        out_sam = "{}M{}I".format(len(seq), nmer)
    else:
        out_sam = "{}M{}I{}M".format(pos, nmer, len(seq) - pos)

    # MD Tag is constant
    out_md = "100"

    # pad reference sequence
    out_ref = ''.join(flatten(
        (itertools.islice(seq, 0, pos), itertools.repeat('_', nmer), itertools.islice(seq, pos, None))))

    # combine everything together and repeat elements as needed!
    return(list(zip(
        out_names,
        itertools.repeat(out_sam, len(nmer_list)),
        itertools.repeat(out_md, len(nmer_list)),
        itertools.repeat(out_ref, len(nmer_list)),
        out_seqs)))

#===============================================================================

def del_len_at(pos, length, seq, viz=False):
    """
    Delete a given length at a certain position in the sequence of interest
    """
    # update to 1 based indexing and ensure you are deleting valid length
    idx = pos - 1
    if idx + length > len(seq):
        raise(ValueError("Length of deletion + position must be less than length of sequence!"))

    # pad the output sequence with _'s
    if viz == True:
        out_seq = ''.join(flatten((
            itertools.islice(seq, 0, idx),
            itertools.repeat('_', length),
            itertools.islice(seq, idx + length, None))))
    else:
        out_seq = ''.join(flatten((
            itertools.islice(seq, 0, idx),
            itertools.islice(seq, idx + length, None))))

    # generate name with form Del_(length)_(@pos)
    out_name = "Del_{}_{}".format(seq[idx:idx + length], pos)

    # generate sam and md tags
    if idx == 0:
        out_sam = "{}D{}M".format(length, len(seq) - length)
        out_md = "0^{}{}".format(seq[idx:idx + length], len(seq) - length)
    elif idx + length >= len(seq):
        out_sam = "{}M{}D".format(len(seq) - length, length)
        out_md = "{}^{}".format(len(seq) - length, seq[idx:idx + length])
    else:
        out_sam = "{}M{}D{}M".format(idx, length, len(seq) - idx - length)
        out_md = "{}^{}{}".format(idx, seq[idx:idx + length], len(seq) - idx - length)

    return((out_name, out_sam, out_md, seq, out_seq,))

#===============================================================================

if __name__ == '__main__':
    # parse input for reference sequence
    parser = argparse.ArgumentParser(
            description='Generate various indels in the specified sequence as determined in the program')
    parser.add_argument('ref', type=str,
            help='sequence to modify (defaults to pasted in)')
    parser.add_argument('-f', '--fasta', dest='fasta', action='store_true',
            help='path to a *.fasta file of the reference sequence')
    parser.add_argument('-v', '--viz', dest='viz', action='store_true',
            help='pad deletions with _ to visuallize better')
    args = parser.parse_args()

    # read a fasta file if req'd
    if args.fasta == True:
        with open(args.ref, 'r') as f:
            next(f)
            ref = f.readline().rstrip()
    else:
        ref = args.ref

    #ref = 'GCTGCCGATTTCCATAAGATGCCTCCACGTCTCCGAAGAACTACATGGTGAATGTGTGAAGGCATTTTGAACCAATCCTCGAGCAGTGTTGTATATATCG'

    # we to insert all 3-mers up to first and last 3 positions
    ins = flatten((insert_nmers_at(x, y, ref)
        for x in flatten((range(0, 4), range(97, 101)))
        for y in range(1, 4)))

    # delete upto the first 5 bases for the first and last 5 positions
    # build a nice iterator to ensure ordering
    order = flatten(zip(itertools.repeat(x, 102 - x), range(1, 102 -x)) for x in range(95, 101))
    end_dels = (del_len_at(x, y, ref, args.viz) for x, y in order)

    front_dels = (del_len_at(x, y, ref, args.viz) for x in range(1, 6) for y in range(1, 6))

    # write to stdout
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['Name', 'CIGAR', 'MD', 'Ref', 'Seq'])
    writer.writerows(list(flatten((ins, front_dels, end_dels))))
