#!/data/home/nlubock/miniconda2/bin/python

"""
bbwrap.py
Nathan Lubock

Simple wrapper for the BBMap suite
"""
import os
import argparse
import subprocess
import multiprocessing as mp
import itertools as it
from functools import partial


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks. Taken from itertools docs
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    """
    args = [iter(iterable)] * n
    return it.izip_longest(*args, fillvalue=fillvalue)

#===============================================================================

def bbtrim(file_pair, ref, out_dir):
    """
    Wrapper function for trimming reads. Note, bbduk.sh uses 7 threads by
    default in the first stage. Should take this into consideration when
    distributing procs.
    Input:
        file_pair - list of paths to the two input reads. Assumes members are
                    of form /path/to/file/my_file_R1.fast(q/a)
        ref - path to ref sequence
    """
    # relevant variables for tuning BBDuk
    k=21 # kmer length for finding contaminants
    mink=8 # kmer length for searching at end of reads
    hdist=2 # hamming distance within kmer match
    hdist2=1 # hamming distance for mink (lower since less number of kmers)
    ktrim='r' # trim reads to the right
    tpe='t' # trim both reads to minimum length of the other
    tbo='f' # trim by overlapping reads
    ram='1g' # constrain memory usage. See bbduk guide for more info
    overwrite='t' # allow overwriting of output files
    threads=5 # number of threads to use after the initial 7 thread stage

    # trim names to current path
    split_names = [os.path.basename(x).split('.') for x in file_pair]
    out_names = [out_dir + x[0] + '.trim.' + x[-1] for x in split_names]

    # NOTE - ASSUMES FILES END IN _R1
    err_base = split_names[0][0][:-3]

    # create string for the os to call
    cmd = 'bbduk.sh in1={} in2={} ref={} k={} mink={} hdist={} hdist2={} ktrim={} overwrite={} tpe={} tbo={} t={} -Xmx{} out1={} out2={}'.format(
            file_pair[0], file_pair[1], ref,  k, mink, hdist, hdist2, ktrim,
            overwrite, tpe, tbo, threads, ram, out_names[0], out_names[1])

    with open(out_dir + err_base + '.trim.err', 'w') as e:
        subprocess.call(cmd, stderr=e, shell=True)
    print 'Finished Trimming: {}!'.format(list(os.path.basename(x) for x in file_pair))

#===============================================================================

def bbfilter(file_pair, ref, out_dir):
    """
    Wrapper function for filtering reads. Note, bbduk.sh uses 7 threads by
    default in the first stage. Should take this into consideration when
    distributing procs.
    Input:
        file_pair - list of paths to the two input reads. Assumes members are
                    of form /path/to/file/my_file_R1.fast(q/a)
        ref - path to ref sequence
    """
    # relevant variables for tuning BBDuk
    k=27 # kmer length for finding contaminants
    hdist=1 # hamming distance within kmer match
    ktrim='f' # delete reads that match
    ram='12g' # constrain memory usage. See bbduk guide for more info
    overwrite='t' # allow overwriting of output files
    threads=1 # number of threads to use after the initial 7 thread stage
    num_n=0 # number of n's to allow in the reads

    # trim names to current path
    split_names = [os.path.basename(x).split('.') for x in file_pair]
    out_names = [out_dir + x[0] + '.filter.' + x[-1] for x in split_names]

    # NOTE - ASSUMES FILES END IN _R1
    err_base = split_names[0][0][:-3]

    # create string for the os to call
    cmd = 'bbduk.sh in1={} in2={} ref={} k={} hdist={} ktrim={} overwrite={} maxns={} t={} -Xmx{} out1={} out2={} stats={}'.format(
            file_pair[0], file_pair[1], ref, k, hdist, ktrim, overwrite, num_n,
            threads, ram, out_names[0], out_names[1],
            out_dir + err_base + '.filter.stats.txt')

    with open(out_dir + err_base + '.filter.err', 'w') as e:
        subprocess.call(cmd, stderr=e, shell=True)
    print 'Finished Filtering: {}!'.format(err_base)

#===============================================================================

def bbqual(merged, level, out_dir):
    """
    Wraper function for filtering merged reads based on their quality scores.

    Parameters:
    -----------
    merged :: [str]
        List of files to filter by quality. Reads should already be merged.
    level :: int
        Reads with average phred values lower than level will be filtered

    Returns:
    --------
    Nothing!

    Output:
    -------
    Fastq file of name $input.phred$level.fastq

    """
    # relevant variables for tuning BBDuk
    overwrite='t' # allow overwriting of output files
    threads=5 # number of threads to use

    # trim input path for output into current directory
    base = os.path.basename(merged).split('.')[0]
    out_base = out_dir + base + '.phred' + str(level)

    cmd = 'bbduk.sh in={} overwrite={} t={} maq={} out={}.fastq'.format(
            merged, overwrite, threads, level, out_base)
    with open(out_dir + out_base + '.err', 'w') as e:
        subprocess.call(cmd, stderr=e, shell=True)
    print 'Finished Filtering: {}!'.format(base)

#===============================================================================

def bbmerge(file_pair, mode, out_dir):
    """
    Wrapper function for merging reads reads.
    Input:
        file_pair - list of paths to the two input reads. Assumes members are
                    of form /path/to/file/my_file_R1.fast(q/a)
        mode - modulate the strictness of the merging
    """
    # relevant variables for tuning BBMerge
    overwrite='t' # allow overwriting of output files
    threads=5 # number of threads to use

    # trim names to current path and remove the _R(1/2)
    # NOTE - assumes that file extensions are the same for both files
    # bbmerge should raise and error if they're not
    split_names = [os.path.basename(x).split('.') for x in file_pair]
    err_base = split_names[0][0][:-3]
    out_name = out_dir + err_base + '.merge.' + split_names[0][-1]

    # create string for the os to call
    if mode != 'perfect':
        cmd = 'bbmerge.sh in1={} in2={} {}=t overwrite={} t={} outm={}'.format(
                file_pair[0], file_pair[1], mode, overwrite, threads, out_name)
    else:
        cmd = 'bbmerge.sh in1={} in2={} pfilter=1 overwrite={} t={} outm={}'.format(
                file_pair[0], file_pair[1], overwrite, threads, out_name)

    with open(out_dir + err_base + '.merge.err', 'w') as e:
        subprocess.call(cmd, stderr=e, shell=True)
    print 'Finished Merging: {}!'.format(err_base)

#===============================================================================

def bbmap(in_file, ref, k, sam, out_dir):
    """
    Wrapper function for mapping reads with bbmap.
    Input:
        in_file - file to map
        ref - path to ref sequence (will cause bbmap to build index in mem)
        k - length of kmer's to use for search (must be same as index if built
            separately!)
        sam - version of sam to output (1.3 vs 1.4)
    """
    # relevant variables for tuning BBMap
    vslow='t' # macro to increase the sensitivity to max
    threads=5 # number of threads to use
    mdtag='t' # compute the mdtag
    ram='8g' # constrain memory usage. See bbduk guide for more info
    maxindel=100 # decrease the maximum length of indels from default 16000

    # handle file names for output
    # NOTE - assumes the only relevant '.' in name is for extensions
    out_base = out_dir + os.path.basename(in_file).split('.')[0]

    # create string for the os to call
    if ref is not None:
        cmd = 'bbmap.sh ref={} in={} k={} vslow={} t={} mdtag={} sam={} maxindel={} -Xmx{} outm={} nodisk'.format(
                ref, in_file, k, vslow, threads, mdtag, sam, ram, maxindel,
                out_base + '.sam')
    else:
        cmd = 'bbmap.sh in={} k={} vslow={} t={} mdtag={} sam={} maxindel={} -Xmx{} outm={}'.format(
                in_file, k, vslow, threads, mdtag, sam, maxindel, ram,
                out_base + '.sam')

    with open(out_base + '.err', 'w') as e:
        subprocess.call(cmd, stderr=e, shell=True)
    print 'Finished Mapping: {}!'.format(out_base)

#===============================================================================

if __name__ == '__main__':
    # create an argument parser a nicer CLI
    parser = argparse.ArgumentParser(
            description='Trim adapter sequences from, remove contaminatns from, merge, or map a list of fast(q/a) files. Note that these files must be paired and should end in _R1 or _R2. Also, filtering contaminants will consume about 10GB per process when searching against the E. Coli genome. NOTE - assumes bbmap.sh and bbduk.sh are in your path!')
    parser.add_argument('-p', '--proc', dest='proc', type=int, default=1,
            metavar='N', choices=range(1, mp.cpu_count()),
            help='number of processors (default=1, max={}). Note that for each processor you specify here, BBDuk will spwan 7 child threads. For example, -p 1 will spawn 7 threads and -p 2 will spawn 14'.format(mp.cpu_count()))

    # add subparser to behave like samtools, e.g. samtools (action) (flags)
    # store subcommand name in action
    subparser = parser.add_subparsers(help='Action to choose', dest='action')
    trim_parse = subparser.add_parser("trim")
    filter_parse = subparser.add_parser("filter")
    merge_parse = subparser.add_parser("merge")
    map_parse = subparser.add_parser("map")
    qual_parse = subparser.add_parser("quality")

    # NOTE:
    # In order to get argparse to handle bash wild cards
    # (e.g. bbwrap.py filter *.fastq) we have to catch all inputs with
    # nargs='+' since Bash will expand wildcards by default. However, we must
    # take care to ensure the positional variables are positioned correctly
    # since nargs='+' will greedily match anything after it

    # trim
    trim_parse.add_argument('ref_path', metavar='adapters',
            help='path to a .fasta file of the adapters to trim')
    trim_parse.add_argument('-i', '--in', dest='file_list', metavar='file', nargs='+',
            help='fast(q/a) file(s) to be trimmed')
    trim_parse.add_argument('-o', '--out', dest='out',
            default='./', type=str,
            help='specify output other than cd')

    # filter
    filter_parse.add_argument('ref_path',
            help='path to the genome(s) of the contaminants to be removed')
    filter_parse.add_argument('-i', '--in', dest='file_list', metavar='file', nargs='+',
            help='fast(q/a) file(s) to be filtered')
    filter_parse.add_argument('-o', '--out', dest='out',
            default='./', type=str,
            help='specify output other than cd')

    # merge
    merge_parse.add_argument('mode',
            default='perfect', type=str,
            choices=['maxloose', 'loose', 'strict', 'maxstrict', 'perfect'],
            help='macros to modulate the strictness of the mering that effectively tunes the false positive rate. "perfect" - only allow perfect overlaps. "maxstrict" - maximally decrease FP rate (allowing for imperfect merges). "maxloose" - maximally increase the merging rate (at the expense of false positives)')
    merge_parse.add_argument('-i', '--in', dest='file_list', metavar='file', nargs='+',
            help='fast(q/a) file(s) to be merged')
    merge_parse.add_argument('-o', '--out', dest='out',
            default='./', type=str,
            help='specify output other than cd')

    # map
    map_parse.add_argument('ref',
            help='path to reference sequence that will be mapped')
    map_parse.add_argument('-i', '--in', dest='file_list', metavar='file', nargs='+',
            help='fast(q/a) file(s) to be mapped')
    map_parse.add_argument('-k', dest='k', default=8, type=int,
            choices=range(8,16), metavar='n',
            help='(Range 8-15)  - kmer length to use during mapping and reference building. Shorter is more sensitive')
    map_parse.add_argument('-S', '--sam', dest='sam', default='1.3', type=str,
            choices=['1.3', '1.4'], metavar='v',
            help='(1.3 or 1.4) - version of sam to output to. Version 1.4 changes the CIGAR string such that mismatches are X and matches are =')
    map_parse.add_argument('-o', '--out', dest='out',
            default='./', type=str,
            help='specify output other than cd')

    # qual
    qual_parse.add_argument('-i', '--in', dest='file_list', metavar='file', nargs='+',
            help='merged fastq file(s) to be filtered')
    qual_parse.add_argument('-c', '--cutoff', dest='cutoff', nargs='?', default=0, metavar='level',
            type=int, choices=range(0,42),
            help='Phred score cut-off (default=0, range=(0,41)). Reads with an average Phred score lower than this value will be discarded.')
    qual_parse.add_argument('-o', '--out', dest='out',
            default='./', type=str,
            help='specify output other than cd')

    #---------------------------------------------------------------------------

    # store everything in args
    args = parser.parse_args()

    # make sure there's an even number of files
    if len(args.file_list) % 2 != 0 and args.action in ['trim', 'filter', 'merge']:
        raise ValueError("Input must have even number of files!")

    # ensure that out paths end in a /
    out_path = os.path.join(args.out, '')

    #---------------------------------------------------------------------------

    # partial the function with the path, then map it over the grouped file_list
    # tune the pool/chunksize for each application
    if args.action == 'trim':
        pool = mp.Pool(args.proc)
        pool.map(partial(bbtrim, ref=args.ref_path, out_dir=out_path),
                grouper(args.file_list, 2), chunksize=1)

    elif args.action == 'filter':
        pool = mp.Pool(args.proc, maxtasksperchild=4)
        pool.map(partial(bbfilter, ref=args.ref_path, out_dir=out_path),
                grouper(args.file_list, 2), chunksize=1)

    elif args.action == 'merge':
        pool = mp.Pool(args.proc)
        pool.map(partial(bbmerge, mode=args.mode, out_dir=out_path),
                grouper(args.file_list, 2), chunksize=1)

    elif args.action == 'map':
        pool = mp.Pool(args.proc)
        pool.map(partial(bbmap, ref=args.ref, k=args.k, sam=args.sam,
            out_dir=out_path), args.file_list)

    elif args.action == 'quality':
        pool = mp.Pool(args.proc)
        pool.map(partial(bbqual, level=args.cutoff, out_dir=out_path),
                args.file_list)
