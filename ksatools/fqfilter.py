#!/usr/bin/env python

import sys
import os
from argparse import ArgumentParser
from Bio import SeqIO
try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

def lesserkmer(s):
    '''returns the lesser of a kmer and its reverse complement'''
    t = revc(s)
    if t < s:
        return t
    else:
        return s


def revc(s):
    '''returns reverse complement of a sequence'''
    intab = "AaCcGgTt"
    outtab = "TtGgCcAa"
    trantab = maketrans(intab, outtab)
    t = s.translate(trantab)[::-1]
    return t


def gccontent(sq):
    '''returns float gc content of sequence'''
    length = len(sq)
    gc_count = 0
    for ch in sq:
        if ch in ['G', 'C', 'g', 'c']:
            gc_count += 1
    try:
        r = float(gc_count) / length
    except:
        r = 0
    return r


def read_index(filename):
    gian = {}
    sys.stderr.write("Processing table %s  ...\n" % (filename,))
    in_idx = open(filename)
    for l in in_idx:
        if l[0] != "#":
            s = l.rstrip().split()
            gian[s[0]] = int(s[1])
    return gian


def kmerabundance(seq, index):
    '''looks up kmer abundance of each kmer in sequence, returns summary statistics'''
    a = []
    for i in range(0, len(seq) - k):
        word = seq[i:i + k]
        w = lesserkmer(word)
        try:
            a.append(index[w])
        except KeyError:
            a.append(0)
    a.sort()
    try:
        median = a[len(a) / 2]
    except IndexError:
        median = 0
    try:
        minimum = a[0]
    except IndexError:
        minimum = 0
    try:
        maximum = a[-1]
    except IndexError:
        maximum = 0
    try:
        average = float(sum(a)) / len(a)
    except IndexError:
        average = 0
    except ZeroDivisionError:
        average = 0
    return (minimum, median, maximum, average)


if __name__ == '__main__':
    usage = "usage: %prog -1 <file1> [-2 <file2>] -i <index>  [-o <outstem>] -l <cutoff>\n  Note: generates outstem.hi.fastq and outstem.lo.fastq"
    parser = ArgumentParser(usage)
    parser.add_argument("-1", "--one", dest="one",
                      required=True, help="Input file 1")
    parser.add_argument("-2", "--two", dest="two",
                      default=None, help="Input file 2 (interleaved if absent)")
    parser.add_argument("-i", "--index", dest="index",
                      required=True, help="input index ")
    parser.add_argument("-t", "--type", dest="typ",
                      default="fastq", help="input datatype (fastq, fasta)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                      default=True, help="Verbose [default off]")
    parser.add_argument("-l", "--cutoff", dest="cutoff",
                      default=None, help="median min coverage cutoff")
    parser.add_argument("-o", "--outstem", dest="outstem",
                      default=None, help="output file stem")
    args = parser.parse_args()
    typ = args.typ
    if not args.typ and args.one[-1] == "a" or args.one[-1] == "A":
        typ = "fasta"
    if not args.one:
        parser.error("Missing input filename")
    if not os.path.isfile(args.one):
        parser.error("Missing input file %s" % args.one)
    if args.two and os.path.isfile(args.two):
        in_two = open(args.two)
    if not args.cutoff:
        sys.stderr.write("Warning: missing cutoff paramter -l\n")
        args.cutoff = 0
    if args.outstem == None:
        args.outstem = args.one
    in_one = open(args.one)
    in_idx = open(args.index)
    if args.verbose or 1:
        sys.stderr.write(
            "Processing sequences %s and table %s  ...\n" % (args.one, args.index))
    sys.stderr.write("Opening output files %s.hi.%s and %s.hi.%s\n" %
                     (args.outstem, typ, args.outstem, typ))
    out_high = open("%s.hi.%s" % (args.outstem, typ), "w")
    out_low1 = open("%s.lo.%s" % (args.outstem, typ), "w")
    giant = {}
    sys.stderr.write("Reading index...\n")
    indexlist = args.index.split(",")
    indexes = []
    for i in range(len(indexlist)):
        giant = read_index(indexlist[i])
        indexes.append(giant)
    k = len(next(indexes[0].iterkeys()))
    sys.stderr.write("Done slurping... set k = %d\n" % (k,))
#   Setup paired-read input
    sys.stderr.write("Looping data: \n")
    records1 = SeqIO.parse(in_one, typ)
    if args.two:
        records2 = SeqIO.parse(in_two, typ)
    else:
        records2 = records1
    n = 0

    for seq_record1 in records1:
        n += 1
        seq_record2 = records2.next()
        seq1 = str(seq_record1.seq)
        seq2 = str(seq_record2.seq)

        if seq1.find("N") == -1 and seq2.find("N") == -1:
            (min1, med1, max1, avg1) = kmerabundance(seq1, indexes[0])
            (min2, med2, max2, avg2) = kmerabundance(seq2, indexes[0])
            seq_record1.description = "%s\tmed%dmer=%d\tmax%dmer=%d\tmin%dmer=%d" % (
                seq_record1.description, k, med1, k, max1, k, min1)
            seq_record2.description = "%s\tmed%dmer=%d\tmax%dmer=%d\tmin%dmer=%d" % (
                seq_record2.description, k, med2, k, max2, k, min2)
            if med1 > float(args.cutoff) and med2 > float(args.cutoff):
                SeqIO.write([seq_record1, seq_record2], out_high, typ)
            else:
                SeqIO.write([seq_record1, seq_record2], out_low1, typ)

    out_low1.close()
    out_high.close()
    if args.verbose:
        sys.stderr.write("Done. \n")
