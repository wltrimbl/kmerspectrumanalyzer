#!/usr/bin/env python 
'''Tool to generate graphs of kmer spectra'''

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import time

def getcolor(a):
    colorlist = ["b", "g", "r", "c", "y", "m", "k"] 
    l = a % len(colorlist)
    return(colorlist[l])


def calcmedian(yd, y, num):
    '''interpolates, returning value of yd corresponding to to num on y'''
    try :
       top = np.max(np.nonzero(y > num)) 
    except:
       top=None
    try :
       bottom = np.min(np.nonzero(y <= num)) 
    except:
       bottom=None
    if top != None and bottom != None:
       cutoff = yd[bottom] + ( ( num - y[bottom] ) / (y[top] - y[bottom]) * (yd[top] -yd[bottom])  )
    elif top != None and bottom == None:
       cutoff = ( ( num *1.0 ) / (y[top] ) * (yd[top])  )
    elif top == None and bottom != None:
       cutoff = yd[bottom]
    else:  
       cutoff = 0
    if num <=1: 
      cutoff=np.ceil(cutoff-.001)
    return(cutoff)

def cleanlabel(label):
    '''Sanitizes graph labels of unintersting file extensions'''
    suffixes=[".histhist", ".fastq", "txt", ".csv", ".037.kmerhistogram"]
    for suffix in suffixes:
        if label.find(suffix) > 0:
            label = label[0:(label.find(suffix))]
    return label

def getmgrkmerspectrum(accessionnumber):
    '''Retrieve kmer spectrum from MG-RAST'''
    import urllib
    import json
    assert accessionnumber[0:3] == "mgm", sys.exit("Data error: field %s not in mgm......... accession number format"%accessionnumber)
    some_url = "http://api.metagenomics.anl.gov/api2.cgi/metagenome_statistics/%s?verbosity=full" % accessionnumber
    if key != None:
        some_url = some_url+"&auth=%s" % key
    sys.stderr.write("Sending request for "+some_url+"\n")
    time.sleep(1)
    try: 
        jsonobject = urllib.urlopen(some_url).read()
    except: 
        sys.stderr.write("Error retrieving %s"%some_url) 
    j = json.loads(jsonobject)
    try: 
        sys.stderr.write("Error with %s\nERROR : %s\n"%(some_url, j["ERROR"]))
        dataarray = None
    except KeyError:
        spectrum = j["qc"]["kmer"]["15_mer"]["data"]
        dataarray = np.array(spectrum, dtype="float")
        dataarray = dataarray[:, 0:2]
    return dataarray

def makegraphs(a, filename, option=6, label=None, n=0):
    '''Draw graphs, one at a time, and add them to the current plot'''
    (cn, c1, yd, yo, zd, zo, y) = calccumsum(a)
    if label == None:
        tracelabel = cleanlabel(filename)
    else: 
        tracelabel = cleanlabel(label)
    b = np.flipud(np.sort(a.view('float,float'), order=['f0'], axis=0).view(np.float))  # sorted by coverage
    c = np.flipud(np.sort(a.view('float,float'), order=['f1'], axis=0).view(np.float))  # sorted by size
    (b_cn, b_c1, b_yd, b_yo, b_zd, b_zo, b_y) = calccumsum(b)
    (c_cn, c_c1, c_yd, c_yo, c_zd, c_zo, c_y) = calccumsum(c)  
    x = np.arange(len(b[:, 0]))                  # rank
    color = getcolor(n)
    if option == 0:
        pA = plt.loglog(cn, c1, "-",  color = color, label=tracelabel)
        pA = plt.loglog(cn, c1, ".",  color = color)
        plt.xlabel("kmer coverage")
        plt.ylabel("number of kmers")
        plt.legend(loc="upper right")
    if option == 0 or option == -1:
        if opts.dump:
            c = np.hstack((cn.reshape((len(cn), 1)), (c1.reshape((len(cn), 1))  )  ))
            sys.stderr.write("saving output table in %s.0.plot.csv\n" % filename) 
            np.savetxt("%s.0.plot.csv" % filename,  c, fmt = ['%d', '%d'] , delimiter="\t" )
    elif option == 1:
        pA = plt.loglog(cn, cn*c1, "-", color=color, label=tracelabel )
        pA = plt.loglog(cn, cn*c1, '.', color=color)
        plt.xlabel("kmer coverage")
        plt.ylabel("kmers observed")
        plt.legend(loc="upper right")
        plt.grid(1)
        if opts.dump:
            c = np.hstack((cn.reshape((len(cn), 1)), ((cn*c1).reshape((len(cn), 1))  )  ))
            sys.stderr.write("saving output table in %s.1.plot.csv\n" % filename) 
            np.savetxt("%s.1.plot.csv" % filename,  c, fmt = ['%d', '%d'] , delimiter="\t" )
    elif option == 2:
        pA = plt.loglog(yo, cn, label=tracelabel )
        plt.xlabel("observed kmers ")
        plt.ylabel("kmer coverage")
        plt.legend(loc="lower left")
        plt.grid(1)
    elif option == 3: 
        y = yo/ yo.max() 
        pA = plt.semilogy(y, cn, label=tracelabel )
        pA = plt.semilogy(y, cn, '.')
        plt.xlabel("fraction of observed sequence")
        plt.ylabel("kmer coverage ")
        plt.grid(1)
        plt.legend(loc="lower left")
    elif option == 4:     # Fraction of distinct kmers vs coverage  NOT RECOMMENDED
        y = yd/ yd.max() 
        pA = plt.semilogy(y, cn, label=tracelabel )
        pA = plt.semilogy(y, cn, '.')
        plt.xlabel("fraction of distinct sequence")
        plt.ylabel("kmer coverage")
        plt.legend(loc="upper right")
        plt.grid(1)
    elif option == 5: 
        z = zo/ zo.max() 
        pA = plt.semilogx( yd, z , '-', color=color )
        pA = plt.semilogx( yd, z, '.', color=color, label=tracelabel )
        plt.xlabel("kmer rank")
        plt.ylabel("fraction of observed kmers")
        plt.xlim((1, 10**9))
        plt.ylim(0, 1)
        plt.grid(1)
        plt.legend(loc="lower left")
    elif option == 6:
        z = zo/ zo.max() 
        pA = plt.loglog( yd, cn, '-', color=color, label=tracelabel)
        pA = plt.loglog( yd, cn, '.', color=color )
        plt.xlabel("kmer rank")
        plt.ylabel("kmer coverage ")
        plt.xlim((1, 10**8))
        plt.ylim(1, 10**7)
        plt.grid(1)
        plt.legend(loc="lower left")
        if opts.dump:
            c = np.hstack((yd.reshape((len(yd), 1)), cn.reshape((len(cn), 1))  )  )
            sys.stderr.write("saving output table in %s.6.plot.csv\n" % filename) 
            np.savetxt("%s.6.plot.csv" % filename,  c, fmt = ['%d', '%d'] , delimiter="\t" )
    elif option == 7:
        pA = plt.plot( x, c_zd, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank")
        plt.ylabel("cuml contig size (bp) ")
        plt.grid(1)
        plt.legend(loc="upper right") 
    elif option == 8:
        pA = plt.plot( x, c_zo/max(c_zo) , '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        plt.legend(loc="upper right") 
    elif option == 9:
        pA = plt.plot( x, c_zo/max(c_zo) , '-', color=color, label=tracelabel)
        plt.xlabel("contig size rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 10:
        pA = plt.plot( x, b_zo/max(b_zo) , '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 11:
        pA = plt.plot( x, b_zo , '.-', color=color,  label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 12:
        pA = plt.plot( c_zd , c_zo , '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")

def calccumsum(a):
    '''Calcaulates the cumulative-sum vectors from a 2d numpy array of [cov, num].  Note depends on upstream sort '''  
    cn = a[:, 0]                          #   Coverage
    c1 = a[:, 1]                          #   number of distinct kmers.
    cp = cn * c1  # elementwise multiply     observed kmers by coverage
    yd = np.flipud(np.flipud(c1).cumsum()) # cumulative number of distinct kmers (top to bottom)
    yo = np.flipud(np.flipud(cp).cumsum()) # cumulative number of observed kmers (top to bottom)
    zd = np.cumsum(c1)                     # cumulative number of distinct kmers (bottom to top)
    zo = np.cumsum(cp)                     # cumulative number of distinct kmers (bottom to top)
    y = zo/ zo.max() 
    return(cn, c1, yd, yo, zd, zo, y)

def printstats(a, filename, filehandle=None, n=0):
    '''Prints summary statistics to filename'''
    cn, c1, yd, yo, zd, zo, y = calccumsum(a)
    T  = zo.max()
    y  = zo/ T
    j  = cn / T
    H = 10**(sum(- c1 * j * np.log(j) /np.log(10)))  # Entropy
    H2 = 1  /(sum( c1*j* j ))                         # Reyni entropy
    w  = yo/yo.max()
    wd = yd
    M90 = calcmedian( wd, w, .9)
    M50 = calcmedian( wd, w, .5)
    M10 = calcmedian( wd, w, .1)
    M100 = calcmedian( wd, w, 1.0)
    F100 = calcmedian( w, wd, 100)
    F10K = calcmedian( w, wd, 10000)
    F1M =  calcmedian( w, wd, 1000000)
    if filehandle == None:
        consensusfh = open(filename, "w")
    else :
        consensusfh = filehandle
    if filehandle == None or n == 0 :
        consensusfh.write( "#filename\tM10\tM50\tM90\tM100\tF100\tF10K\tF1M\tH\tH2\n")
    consensusfh.write( "%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%.1f\t%.1f\n" % 
                      (filename, M10, M50, M90, M100, F100, F10K, F1M, H, H2))
    if filehandle == None:
        consensusfh.close()
 
def loadfile(filename):
    '''Loads file'''
    try: 
        if filename.find("stats.txt") >=0:  # velvet contig stats format
            a = np.loadtxt(filename, usecols=(5, 1), skiprows=1)
        else:
            a = np.loadtxt(filename)        # default bare-bones spectrum format 
        return a 
    except IOError:
        sys.stderr.write("ERROR: Can't find file %s\n" % filename)
        return None

def main(filename, option=6, label=None, n=0 ):
    '''Main driver; loads file and invokes makegraphs and printstats 
    to append graphics from each file onto the figure'''
    logfh = open(opts.logfile, "a")
    if opts.filetype.upper() == "MGM":
        a = getmgrkmerspectrum(filename)
    elif opts.filetype == "file":
        a = np.atleast_2d(loadfile(filename))
    else: 
        raise ValueError("%s is invalid type (valid types are mgm and file)"%opts.filetype ) 
    if label == None:
        label = filename
    if a != None:
        a = (a[np.lexsort((a[:, 1], a[:, 0]))])
        sys.stderr.write("Making graphs for %s\n" % filename)
        makegraphs(a, filename, option, label, n=n )
        try: 
            sys.stderr.write("Printing stats in logfile %s %d\n" % (opts.logfile, n))
            printstats(a, filename, filehandle=logfh, n=n )
            printstats(a, filename, filehandle=sys.stdout, n=n) 
            n += 1
        except: 
            sys.stderr.write("Error printing stats for %s\n" % filename)
            print "Unexpected error:", sys.exc_info()[0]
    else:
        sys.stderr.write("Error with dataset %s\n" % filename)
    return n  
if __name__ == '__main__':
    usage  = "usage: %prog [options] <datafile> [<datafile2> <datafile3...]"
    parser = OptionParser(usage)
    parser.add_option("-d", "--dump",   dest="dump",   action="store_true", 
         default=False, help="dump table with outputs ")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", 
         default=False, help="verbose")
    parser.add_option("-o", "--outfile", dest="outfile", action="store", 
         default=None, help="dump table with outputs ")
    parser.add_option("-g", "--graph",  dest="option", action="store", type="int", 
         default="6", help="Graph number ")
    parser.add_option("-i", "--interactive", dest="interactive", action="store_true", 
         default=False, help="interactive mode--draw window")
    parser.add_option("-l", "--list", dest="filelist",  
         default=None, help="file containing list of targets and labels")
    parser.add_option("-t", "--type", dest="filetype", 
         default="file", help="type for file list (file,mgm)")
    parser.add_option("-a", "--appendlogfile", dest="logfile", 
         default="kmers.log", help="logfile for summary statistics")
  
    (opts, args) = parser.parse_args()
    option = opts.option
    if opts.outfile: 
        imagefilename = "%s.%d.pdf" % (opts.outfile, option)
    else: 
        imagefilename = "out.%d.pdf" % (option,)
        sys.stderr.write("Warning, using default filename %s\n" % (imagefilename,))
    if opts.filetype == "mgm":  
        try:
            key = os.environ["MGRKEY"]
        except KeyError:
            key = ""

    graphcount = 0
    if opts.filelist: 
        assert os.path.isfile(opts.filelist), "File %s does not exist"%opts.filelist
        IN_FILE = open(opts.filelist, "r") 
        for line in IN_FILE:
            a = line.rstrip().split("\t")
            if len(a) == 1: 
                a.append(a[0])
            sys.stderr.write( "%s  %s \n" % (a[0], a[1]) ) 
            graphcount = main(a[0], option, label=a[1], n=graphcount)
    else:
        for f in args :
            filen = f
            graphcount = main(filen, option, n=graphcount)
    if option != -1:
        sys.stderr.write("Writing graph into file %s\n" % (imagefilename))
        plt.savefig(imagefilename)
    if opts.interactive:
        plt.show()
    else:
        sys.stderr.write( "Use -i to open widow with graph\n")
    if graphcount == 0 :
        sys.stderr.write("ERROR:  no data found!\n") 