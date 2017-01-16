# bowtieWrap.py
# Nathan Lubock
#
#
# USAGE:
#      python bowtieWrap.py /path/to/file /path/to/bt2_index

import os
import sys
import subprocess
import argparse
import multiprocessing as mp
from functools import partial


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks. Taken from itertools docs
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    """
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)

#=============================================================================

def bowtie(in_file, ref_path):
    """
    Generate a bowtie cmd with some standard options to debug
    Input:
    Output:
        cmd - bash command to run bowtie2
    """
    # default variables
    threads = 5
    seedLen = 5 # length of seed substring during multiseed alignment
    MM = 1 # num(mismatches) during multiseed alignment stage
    interval = 'C,1' # interval between seeds during multiseed alginment
    scoreMin = 'C,-200' # min score to report
    gaps = 1 # disallow gaps within the _ bp of a read
    rdg = '5,1' # read gap open penalty A + N*B (default 5,3)

    # -x = basename for reference genome
    # -U = name for unpaired (or merged) fastq's
    # --un = reads that failed to align
    # -S = sam to write to (instead of stdout)

    # handle file names for output
    # NOTE - assumes the only relevant '.' in name is for extensions
    out_base = os.path.basename(in_file).split('.')[0]

    # generate your command!
    cmd = ("bowtie2 --end-to-end --time --no-unal -L {} -N {} -i {} --gbar {} --score-min {} -p {:g} --rdg {} -x {} -U {} --un {} -S {}").format(
           seedLen, MM, interval, gaps, scoreMin, threads, rdg, ref_path,
           in_file, out_base + '.unalign.fastq', out_base + '.sam')
    with open(out_base + '.err', 'w') as e:
        subprocess.call(cmd, stderr=e, shell=True)
    print('Finished Mapping: {}'.format(out_base))


#=============================================================================


# Main body
# Usage:
#       bowtieWrap.py $file_list $bt2_index
# where:
#       $file_list = file containing paths tho paired end reads
#       $bt2_index = full path specifying basename of bt2 index
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description='Wrapper for distributing lots of bowtie2 jobs')
    parser.add_argument('-p', '--proc', dest='proc', type=int, default=1,
            metavar='N', help='number of processors (default=1)')
    parser.add_argument('file_list',
            help = 'a plain-text file specifying the path to the fastq file(s)')
    parser.add_argument('ref_path',
            help = 'path to bowtie2-index. Do not provide an extension')
    args = parser.parse_args()

    with open(args.file_list, 'r') as f: # inefficient for large lists
        file_list = f.read().splitlines()

    # Set up multiprocessing
    # Be sure to adjust the number of threads accordingly if you add workers
    pool = mp.Pool(args.proc)
    pool.map(partial(bowtie, ref_path=args.ref_path), file_list)

