#!/usr/bin/env python2

# filter_fastq.py

# filter sequence ids by identifiers in a file.

# SNF March 2016

# Changelog
#
# (2021-12-11)
# 1. Modified script at https://github.com/Floor-Lab/filter-fastq to provide keywords (in regex format) to filter fastq sequences.
# 2. Tweaks to speed up the process.

import sys
import argparse
import re
import gzip
import os
import threading
import time
from itertools import islice, izip
import multiprocessing

parser = argparse.ArgumentParser(description="Filter fastq file")

parser.add_argument(
    "-i",
    "--input",
    help="File to filter unpaired reads from (default stdin)",
    default="stdin",
)
parser.add_argument("-1", "--read1", help="Paired-end fastq file one")
parser.add_argument("-2", "--read2", help="Paired-end fastq file two")
parser.add_argument(
    "-p",
    "--num-threads",
    help="Number of threads to start (currently for PE reads only)",
    type=int,
)
parser.add_argument(
    "-o",
    "--output",
    help="Output filename (default stdout; interpreted as a basename for PE reads)",
    default="stdout",
)
parser.add_argument(
    "-f",
    "--filter_file",
    help="File containing strings to include based on IDs in fastq",
)
parser.add_argument(
    "-k",
    "--keyword_search",
    help="FILTER_FILE contains keywords (as regex) to filter input fastq file(s)",
    action="store_true",
)
parser.add_argument(
    "-v",
    "--invert",
    help="Invert match (exclude matching entries)",
    action="store_true",
)
parser.add_argument(
    "--gzip", help="gzip compress the output", action="store_true")

args = parser.parse_args()

if args.input != "stdin" and (args.read1 or args.read2):
    sys.exit(
        "Fatal: please provide only one of --input or --read1/--read2 but not both"
    )

if (args.read1 and not args.read2) or (args.read2 and not args.read1):
    sys.exit(
        "Fatal: both --read1 and --read2 must be provided when processing paired-end reads"
    )

if args.read1:  # reads are paired now if this is true
    paired = True
else:
    paired = False

if not paired:
    if args.output == "stdout" or args.output == "-":
        outfile = sys.stdout
    else:
        if args.gzip:
            outfile = gzip.open(args.output, "wb")
        else:
            outfile = open(args.output, "w")

    if args.input == "stdin" or args.input == "-":
        infile = sys.stdin
    else:
        filename, file_extension = os.path.splitext(args.input)
        if file_extension == ".gz":
            infile = gzip.open(args.input, "rb")
        elif file_extension in [".fastq", ".fq"]:
            infile = open(args.input, "r")
        else:
            sys.exit("Fatal: unknown file extension for filename %s" %
                     args.input)

else:  # paired
    filename, file_extension = os.path.splitext(args.read1)
    if file_extension == ".gz":
        infile_read1 = gzip.open(args.read1, "rb")
    elif file_extension in [".fastq", ".fq"]:
        infile_read1 = open(args.read1, "r")
    else:
        sys.exit("Fatal: unknown file extension for filename %s" % args.read1)

    filename, file_extension = os.path.splitext(args.read2)
    if file_extension == ".gz":
        infile_read2 = gzip.open(args.read2, "rb")
    elif file_extension in [".fastq", ".fq"]:
        infile_read2 = open(args.read2, "r")
    else:
        sys.exit("Fatal: unknown file extension for filename %s" % args.read2)

    if args.gzip:
        outfile_read1 = gzip.open(args.output + ".1.fastq.gz", "wb")
        outfile_read2 = gzip.open(args.output + ".2.fastq.gz", "wb")
    else:
        outfile_read1 = open(args.output + ".1.fastq", "w")
        outfile_read2 = open(args.output + ".2.fastq", "w")

    sys.stderr.write(
        "Filtering reads from paired files %s and %s.\n" % (
            args.read1, args.read2)
    )

filter_list = {i.strip(): re.compile(i.strip())
               for i in open(args.filter_file)}

sys.stderr.write("Read %d identifiers to filter.\n" % len(filter_list))

# note this function uses filter_list as defined globally earlier. can't pass it in b/c need one argument for pool.imap


def filter_pair(lines):
    # lines is a set of four lines from both files that are stitched together. Order is:
    #   file1line1, file2line1, file1line2, file2line2, etc...
    seqheader1 = lines[0]
    seqid1 = seqheader1.split()[0]

    seqheader2 = lines[1]
    seqid2 = seqheader2.split()[0]

    assert seqid1 == seqid2, "Out of sync pairs detected with ids %s %s" % (
        seqid1,
        seqid2,
    )

    if seqheader1 == seqheader2:
        sys.stderr.write("Same sequence header shared by reads %s %s" % (
            seqheader1,
            seqheader2,
        ))

    assert (lines[4] == lines[5]) and (
        lines[4] == '+\n'), "Incorrect input file format"

    return_status = False
    if args.keyword_search:
        for keyword in filter_list:
            if ((filter_list[keyword].search(seqheader1)) or (filter_list[keyword].search(seqheader2))):
                return_status = True
                break
        if args.invert:
            return_status = not return_status
        if return_status:
            return lines
        return return_status
    else:
        if seqid1 in filter_list:
            return_status = True
        if args.invert:
            return_status = not return_status
        if return_status:
            return lines
        return return_status


def filter_read(lines):
    # lines is a set of four lines from the file.
    seqheader = lines[0]
    seqid = seqheader.split()[0]

    return_status = False
    if args.keyword_search:
        for keyword in filter_list:
            if (filter_list[keyword].search(seqheader)):
                return_status = True
                break
        if args.invert:
            return_status = not return_status
        if return_status:
            return lines
        return return_status
    else:
        if seqid in filter_list:
            return_status = True
        if args.invert:
            return_status = not return_status
        if return_status:
            return lines
        return return_status


n_reads = 0
n_filtered = 0
pool = multiprocessing.Pool(processes=args.num_threads)

if paired:
    # the izip below makes an iterator that loops over both files in sets of four. should be parallel-safe
    res = pool.imap(filter_pair, izip(*[infile_read1, infile_read2] * 4))
else:
    # one-file iterator
    res = pool.imap(filter_read, izip(*[infile] * 4))

while 1:
    if not n_reads % 1000000:
        sys.stderr.write("Processed %d reads.\n" % n_reads)
    try:
        lines = res.next()
        n_reads += 1
        if lines:
            if paired:
                outfile_read1.write(lines[0] + lines[2] + lines[4] + lines[6])
                outfile_read2.write(lines[1] + lines[3] + lines[5] + lines[7])
            else:
                outfile.write("".join(lines))
            if not args.invert:
                n_filtered += 1
        else:  # don't output this
            if args.invert:
                n_filtered += 1

    except StopIteration:
        break

pool.close()
pool.join()

sys.stderr.write(
    "Processed %d reads, filtered %d reads. (%.2f%%).\n"
    % (n_reads, n_filtered, n_filtered / float(n_reads))
)
