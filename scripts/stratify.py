#!/usr/bin/env python3
'''Tool to generate bar graphs of data-depth fractions and
basepair-by-depth distributions from kmer spectra
'''

from __future__ import print_function
import sys
import os
import argparse
import numpy as np
import matplotlib as mpl

from ksatools.ksatools import loadfile, stratify, cleanlabel


def plotstratifiedfracs(inlabels, inspectra):
    '''Given list of labels and list of spectra, splits spectra
    up and produces bar plot of the stratified spectra's usage
    fractions
    '''
    colors = ["#CCFFCC", "#99FF99", "#55FF55", "#33CC33", "#009933", "#003300"]
    bands = []
    fracs = []
    sizes = []
    plt.grid(linestyle="-", axis="x", zorder=-10)
    cleanlabels = [cleanlabel(l) for l in inlabels]
    for i in range(len(inlabels)):
        label = cleanlabels[i]
        spectrum = inspectra[i]
        sys.stderr.write("Stratifying " + label + "...\n")
        BANDS = [1, 3, 30, 300, 3000, 30000, 30000000]
        band, frac, size = stratify(spectrum, bands=BANDS)
        bands.append(band)
        fracs.append(frac)
        sizes.append(size)
    for l in range(len(inlabels)):
        for i in range(len(bands[0]) - 1):
            if l == 0:
                plt.barh(l, (fracs[l][i] - fracs[l][i + 1]),
                         left=(fracs[l][i + 1]), color=colors[i],
                         label=str(bands[0][i]) + "-" + str(bands[0][i + 1]),
                         alpha=1.0, zorder=0)
            else:
                plt.barh(l, (fracs[l][i] - fracs[l][i + 1]),
                         left=(fracs[l][i + 1]), color=colors[i],
                         alpha=1.0, zorder=0)
    pos = np.arange(len(inlabels)) + 0.0
    plt.yticks(pos, cleanlabels)
    plt.xlabel("Cumulative data fraction")
    plt.tight_layout()
    if not opts.suppresslegend:
        plt.legend(loc="upper left")
    plt.show()


def plotstratifiedsizes(inlabels, inspectra):
    '''Given list of labels and spectra, produces stacked bar graphs
    on a log scale of the cumulative number of basepairs above or
    equal to each depth boundary.
    '''
    colors = ["#CCFFCC", "#99FF99", "#55FF55", "#33CC33", "#009933", "#003300"]
    plt.grid(linestyle="-", axis="x", zorder=-10)
    bands = []
    fracs = []
    sizes = []
    BANDS = [1, 3, 30, 300, 3000, 30000, 30000000]
    cleanlabels = [cleanlabel(l) for l in inlabels] 
    maxsize = 0
    for i in range(len(inlabels)):
        label = cleanlabels[i]
        spectrum = inspectra[i]
        sys.stderr.write("Stratifying " + label + "...\n")
        band, frac, size = stratify(spectrum, bands=BANDS)
        if np.max(size) > maxsize:
            maxsize = np.max(size)
        bands.append(band)
        fracs.append(frac)
        sizes.append(size)
    for l in range(len(inlabels)):
        for i in range(len(bands[0]) - 1):
            sizec = np.array(sizes[l])
            if l == 0:
                plt.barh(l, (sizec[i + 1] - sizec[i]), left=(sizec[i]),
                         color=colors[i], label=str(bands[0][i]) +
                         "-" + str(bands[0][i + 1]), log=True)
            else:
                plt.barh(l, (sizec[i + 1] - sizec[i]), left=(sizec[i]),
                         color=colors[i], log=True)
    pos = np.arange(len(inlabels)) + 0.0
    xlimit = max(1E9,  10**int(np.log10(maxsize) + 1))
    plt.xlim((1, xlimit))  # lower xlimit should be an even power of 10
    plt.yticks(pos, cleanlabels)
    plt.xlabel("Distinct kmers (basepairs)")
    plt.tight_layout()
    if not opts.suppresslegend:
        plt.legend()
    plt.show()


def plotstratifiedcomplexity(inlabels, inspectra):
    '''Given list of labels and spectra, produces stacked bar graphs
    of cumulative basepairs explained stratified by bands of complexity
    '''
    colors = ["#FFCCCC", "#FF9999", "#FF5555", "#CC3333", "#990033", "#330000","#CC3333", "#990033", "#330000"]
    colors = [ "#330000", "#990033", "#992233", "#CC3333", "#CC5533", "#FF5555", "#FF7755", "#FF9999", "#FFCCCC", "#FFFFCC", ]
    colors = [ "black", "red", "orange", "yellow", "green", "blue", "violet", "gray", "white"]
    colors = [ "#000000", "#220044", "#880000", "#CC0044", "#FF0000", "#FF4488", "#FF8888", "#FFCCFF", "#FFFFFF", "#888888", ]
    plt.grid(linestyle="-", axis="x", zorder=-10)
    bands = []
    fracs = []
    sizes = []
    BANDS = [100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000]
    BANDLABELS = ['100', '1k', '10k', '100k', '1M', '10M', '100M', '1G', '10G']
    cleanlabels = [cleanlabel(l) for l in inlabels]
    for i in range(len(inlabels)):
        label = cleanlabels[i]
        spectrum = inspectra[i]
        sys.stderr.write("Stratifying " + label + "...\n")
        band, frac, size = stratify(spectrum, bands=BANDS, bandtype="complexity")
        bands.append(band)
        fracs.append(frac)
        sizes.append(size)
    for l in range(len(inlabels)):
        for i in range(len(bands[0]) - 1):
            if l == 0:
                plt.barh(l, (fracs[l][i] - fracs[l][i + 1]),
                         left=(fracs[l][i + 1]), color=colors[i],
                         label=str(bands[0][i]) + "-" + str(bands[0][i + 1]),
                         alpha=1.0, zorder=0)
            else:
                plt.barh(l, (fracs[l][i] - fracs[l][i + 1]),
                         left=(fracs[l][i + 1]), color=colors[i],
                         alpha=1.0, zorder=0)
    pos = np.arange(len(inlabels)) + 0.0
    plt.xlim((0, 1))
    plt.yticks(pos, cleanlabels)
    plt.xlabel("Cumulative data fraction")
    plt.tight_layout()
    if not opts.suppresslegend:
        plt.legend()
    plt.show()


def summarizestrata(inlabels, inspectra):
    '''Prints one-line table-style summary of cumulative fractions and
    sizes
    '''
    BANDS = [1, 3, 10, 30, 100, 300, 1000, 3000,
             10000, 30000, 100000, 300000, 1000000]
    for spectrum, label in zip(inspectra, inlabels):
        band, frac, size = stratify(spectrum, bands=BANDS)
        bandsbanner = map(str, BANDS)
        print("#name\t" + "\t".join(bandsbanner) +
              "\t" + "\t".join(bandsbanner))
        print(label + "\t" + "\t".join(map(str, frac)) + "\t",)
        print("" + "\t".join(map(str, size)))
    return


if __name__ == '__main__':
    usage = '''stratify.py [options] <list of target files, labels>
'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true",
        default=False, help="verbose")
    parser.add_argument(
        "-o", "--outfile", dest="outfile", action="store",
        default=None, help="output file name")
    parser.add_argument(
        "-g", "--graph", dest="option", action="store", type=int,
        default="0", help="graph type 0 = fraction, 1 = basepairs, 2 = fraction by basepairs, -1 = table")
    parser.add_argument(
        "-i", "--interactive", dest="interactive", action="store_true",
        default=False, help="interactive mode--draw window")
    parser.add_argument(
        "-l", "--list", dest="filelist",
        default=None, help="file containing list of targets and labels")
    parser.add_argument(
        "-t", "--type", dest="filetype",
        default="file", help="type for file list (file,mgm)")
    parser.add_argument(
        "-w", "--writetype", dest="writetype",
        default="pdf", help="file type for output (pdf,png)")
    parser.add_argument(
        "-s", "--suppresslegend", dest="suppresslegend", action="store_true",
        default=False, help="supress display of legend")
    parser.add_argument(
        "-n", "--name", dest="title",
        default=None, help="Name for graph, graph title")

    opts = parser.parse_args()
    writetype = opts.writetype
    if opts.filelist is None:
        sys.exit(usage)
    filelist = opts.filelist
    if filelist is None:
        sys.exit(usage)
    assert writetype == "png" or writetype == "pdf" or writetype == "eps"
    if opts.outfile:
        imagefilename = opts.outfile
    else:
        imagefilename = opts.filelist + "." + opts.writetype
        sys.stderr.write("Warning, using default filename %s\n" %
                         (imagefilename,))
    # only invoke interactive backend if requested with -i
    # this stabilizes behavior on non-interactive terminals
    if not opts.interactive:
        mpl.use("Agg")
    else:
        mpl.use('TkAgg')
    import matplotlib.pyplot as plt

    graphcount = 0
    # Loop over input identifiers, and skip if main()
    # fails to produce some traces
    spectra = []
    labels = []
    if opts.filelist:
        if not os.path.isfile(opts.filelist):
            sys.exit("File {} does not exist".format(opts.filelist))
        IN_FILE = open(opts.filelist, "r")
        for line in IN_FILE:
            if line[0] != "#":
                a = line.strip().split("\t")
                if len(a[0]) > 0:
                    if len(a) == 1:
                        a.append(a[0])
                    sys.stderr.write("%s\t%s\n" % (a[0], a[1]))
                    contents = loadfile(a[0])
                    if contents != []:
                        spectra.append(contents)
                        labels.append(a[1])
    if opts.option == 0:
        plotstratifiedfracs(labels, spectra)
    elif opts.option == 1:
        plotstratifiedsizes(labels, spectra)
    elif opts.option == 2:
        plotstratifiedcomplexity(labels, spectra)
    elif opts.option == -1:
        summarizestrata(labels, spectra)
        sys.exit()
    else:
        sys.exit("Valid graph types are -g 0, 1, and -1")
    if opts.interactive:
        sys.stderr.write("Interactive mode (no file output)\n")
        plt.show()
    else:
        sys.stderr.write("Use -i to open widow with graph\n")
        sys.stderr.write("Writing graph into file %s\n" % (imagefilename))
        plt.savefig(imagefilename)
