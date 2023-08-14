import pandas as pd
import pysam
import numpy as np
from scipy.stats import binom
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt
import statsmodels.api as sm
import matplotlib.pyplot as plt
import sys
import argparse

class ase_simpler():
    def __init__(self, args):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--inputBamFile", dest="inputBam")
        parser.add_argument("-h", "--hetSites", dest="hetSites")
        self.inputBam = args[0]
        if (self.inputBam[-4:] != '.bam'):
            print("Error Input Bam File")
            exit()
        self.hetSites = args[1]

    # summary: reads in read counts and locations of snps
    # parameter (inputBam): input bam file
    # parameter (hetSites): location for hetSites, given format "chr | start | end"
    # returns: return a pandas df of chr, index, and basecounts
    def readText(inputBam, hetSites):
        locationFile = open(hetSites, "r").read()
        samfile = pysam.AlignmentFile(inputBam, "rb")
        dfDict = {'chr': [], 'start': [], 'end': [], 'A': [], 'C': [], 'G': [], 'T': []}
        for line in locationFile:
            baseDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            # extract information from location file
            lineLst = line.split(" | ")
            if len(lineLst) != 3:
                print("Improper het site location format")
                exit()
            try:
                chr = lineLst[0]
                start = int(lineLst[1])
                end = int(lineLst[2])
            except TypeError:
                print("Not a valid het site index.")
                exit()
            # initialize pileup object for the chr, start, and end extracted
            pileupColumn = samfile.pileup(chr, start, end)
            # iterate through the pileup 
            for pileupRead in pileupColumn.pileups:
                # base at specified position
                base = pileupRead.alignment.query_sequence[pileupRead.query_position]
                # add to the count of bases that appear for that specific site
                baseDict[base] += 1
            dfDict['chr'].append(chr)
            dfDict['start'].append(start)
            dfDict['end'].append(end)
            dfDict['A'].append(baseDict['A'])
            dfDict['C'].append(baseDict['C'])
            dfDict['G'].append(baseDict['G'])
            dfDict['T'].append(baseDict['T'])
        return pd.DataFrame(dfDict)
        
    def run(self, args):

    
aseModel = ase_simpler(sys.argv)
aseModel.run()