import pandas as pd
import pysam
import numpy as np
from scipy.stats import binom
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import sys
import argparse

# written for python 2

class ase_simpler():
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--inputBamFile", dest="inputBam", help="input bam file from alignment")
        parser.add_argument("-hs", "--hetSites", dest="hetSites", help="input het site information, in the format \"chr | start | end\"")
        parser.add_argument("-rc", "--readCounts", dest="readCounts", help="get read counts from the read counts file")
        parser.add_argument("-v", "--vcf", dest="vcf", help="get het site information using the vcf")
        parser.add_argument("-p", "--peaks", dest="peaks", help="get peaks with a start and end index")
        args = parser.parse_args()
        self.inputBam = args.inputBam
        self.hetSites = args.hetSites
        self.readCounts = args.readCounts
        self.vcf = args.vcf
        self.peaks = args.peaks

    # summary: reads in read counts and locations of snps
    # parameter (inputBam): input bam file
    # parameter (hetSites): location for hetSites, given format "chr | index"
    # returns: writes into a csv file
    def readText(self):
        inputBam = self.inputBam
        hetSites = self.hetSites
        locationFile = open(hetSites, "r").readlines()
        samfile = pysam.AlignmentFile(inputBam, "rb")
        outputFile = open("readCounts.csv", "w")
        dfDict = {'chr': [], 'index': [], 'A': [], 'C': [], 'G': [], 'T': []}
        count = 0
        outputFile.write("chr,index,A,C,G,T")
        for line in locationFile:
            if not count % 100:
                print(str(count) + " : " + str(len(locationFile)))
            baseDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            # extract information from location file
            lineLst = line.split(" | ")
            chr = lineLst[0]
            index = int(lineLst[1])
            # initialize pileup object for the chr and index extracted
            for pileupColumn in samfile.pileup(chr, index, index+1):
                # iterate through the pileup
                for pileupRead in pileupColumn.pileups:
                    # base at specified position
                    if (pileupRead.query_position == None):
                        continue
                    base = pileupRead.alignment.query_sequence[pileupRead.query_position]
                    # add to the count of bases that appear for that specific site
                    baseDict[base] += 1
            count += 1
            outputFile.write("{0},{1},{2},{3},{4},{5}\n".format(chr,index,baseDict['A'],baseDict['C'],baseDict['G'],baseDict['T']))
        return pd.DataFrame(dfDict)

    # summary: create a compatible vcf as a pandas dataframe for reference and alternate information
    # parameter (vcf_file): a vcf file with the header
    # returns: a DataFrame with all the VCF information
    def proc_vcf(vcfFile):
        vcf = open(vcfFile, "r")
        vcf_cut = open('{0}_abbrev.txt'.format(vcfFile), 'w')
        for line in vcf:
            if (line[0] == "#"):
                continue
            vcf_cut.write(line + "\n")
        vcf_cut.close()
        vcf_df = pd.read_csv("{0}_abbrev.txt".format(vcfFile), sep='\t')
        return vcf_df
    
    # summary: get a pileup read counts dataframe
    # parameter: read counts file
    # returns: a DataFrame with all the read counts information
    def getReadCounts(self):
        readCounts = self.readCounts
        df = pd.read_csv(readCounts, sep="\t")
        return df

    # summary: get complimentary base
    # parameter: character base pair
    # returns: complimentary base pair
    def comp(char):
        if char == "A":
            return "T"
        elif char == "T":
            return "A"
        elif char == "C":
            return "G"
        elif char == "G":
            return "C"
        else:
            return ""

    # summary: creates a dictionary in the format: {index: (refCount, altCount), ...}
    # parameter: takes in a dataframe of the read counts at each index and dataframe of the vcf
    # returns: the dictionary of the format: {index: (refCount, altCount), ...}
    def ref_alt(df, vcf_df):
        # reference sequence string placeholder
        refseq = ""
        # alternate sequence string placeholder
        altseq = ""
        # alternate count
        alt = 0
        # reference count
        ref = 0
        # length of the vcf (number of snps in the sample)
        n = len(df.index)
        # dictionary to store the information
        ref_alt = {}
        # erroneous samples go in a list (no counts at all)
        err = []
        # ignore indels
        indelCount = 0
        # iterate through all of the cases
        for i in range(n):
            # set refseq and altseq to the items in each if they are an alphanumeric characte
            refseq = "".join(c for c in vcf_df.iloc[[i]]["REF"].item() if c.isalnum())
            altseq = "".join(c for c in vcf_df.iloc[[i]]["ALT"].item() if c.isalnum())
            # indel events actually have strings of longer than 1, therefore just continue and keep track
            if len(refseq) > 1 or len(altseq) > 1:
                indelCount += 1
                continue
            # reference and alternate counts are gathered by using this dictionary referencing method
            ref = df.iloc[[i]][refseq.upper()].item() + df.iloc[[i]][comp(refseq).lower()].item()
            alt = df.iloc[[i]][altseq.upper()].item() + df.iloc[[i]][comp(altseq).lower()].item()
            # if the counts are both 0, then clearly there is an issue (snp with no read counts)
            if ref == 0 and alt == 0:
                err.append([ref, alt])
                continue
            # if the vcf_index is not a key then create a list
            if vcf_df.iloc[[i]]["POS"].item() not in ref_alt.keys():
                ref_alt[vcf_df.iloc[[i]]["POS"].item()] = []
            # append the read counts tuple to the list
            ref_alt[vcf_df.iloc[[i]]["POS"].item()].append((ref, alt))
        return ref_alt, err

    # summary: generates a DataFrame with the start and end of the peaks from the peaks file
    # parameters: file of the peaks formatted: "CHR START END ID SAMPLE1 ..." tab separated, and the delimiter
    # returns: generated DataFrame
    def gen_peaks_df(peaksFile, sep=","):
        peaks = pd.read_csv(peaksFile, sep=sep)
        size = peaks.shape[0]
        d = {}
        for i in range(size):
            d[i] = [peaks.iloc[[i]]["CHR"].item(), str(int(peaks.iloc[[i]]["START"].item() - 150)), str(int(peaks.iloc[[i]]["END"].item() - 150))]
        peaks_df = pd.DataFrame.from_dict(d, orient="index", columns=["CHR", "START", "END"])
        return peaks_df

    # summary:

    # runs the main method
    def run(self):
        # get DataFrame of the read counts piled up at each het site
        pileupDf = self.getReadCounts()
        # get a dataframe containing the vcf information
        vcf_df = self.proc_vcf(self.vcf)
        # get a dictionary of the format
        ref_alt, err = self.ref_alt(pileupDf, vcf_df)
        # get peaks from the peaks file
        peaks = self.gen_peaks_df(self.peaks, "\t")


aseModel = ase_simpler()
aseModel.run()
