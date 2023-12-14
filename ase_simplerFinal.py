import pandas as pd
import pysam
import numpy as np
from scipy.stats import binom
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import sys
import argparse
from pathlib import Path
import timeit
import os

# written for python 2
# python ase_simplerFinal.py -rc G1356B_S28_L004output.txt -v final288.all_chrm.vcf -p quant_peaks_1FDR_minOverlap2_300bpExt_150_TMM_no_rep_outliers.mid_peak.288.txt
# python ase_simplerFinal.py -rc ReadCounts -v allChrmVCF -p peaks

# 10/05/2023
# only look into autosomes
# leave out indels
class ase_simpler():
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--inputBamFile", dest="inputBam", help="input bam file from alignment")
        parser.add_argument("-hs", "--hetSites", dest="hetSites", help="input het site information, in the format \"chr | start | end\"")
        parser.add_argument("-rc", "--readCounts", dest="readCounts", help="get read counts from the read counts file")
        parser.add_argument("-v", "--vcf", dest="vcf", help="get het site information using the vcf")
        parser.add_argument("-p", "--peaks", dest="peaks", help="get peaks with a start and end index")
        parser.add_argument("-o", "--outputDirectory", dest="outDir", help="output directory for pVals dictionary (outputs csv)")
        parser.add_argument("-ch", "--chromosome", dest="chrom", help="chromosome number")
        args = parser.parse_args()
        self.inputBam = args.inputBam
        self.hetSites = args.hetSites
        self.readCounts = args.readCounts
        self.vcf = args.vcf
        self.peaks = args.peaks
        self.outDir = args.outDir
        self.chrom = args.chrom
        isExist = os.path.exists(self.outDir)
        if not isExist:
            os.makedirs(self.outDir)

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
    def proc_vcf(self, vcfFile):
        vcf = open(vcfFile, "r")
        filename = Path(vcfFile).stem
        vcf_cut = open('{0}_abbrev.txt'.format(filename), 'w')
        header = ""
        for line in vcf:
            if (line[0] == "#" and line[1:6] == "CHROM"):
                vcf_cut.write(line[1:] + "\n")
                continue
            elif (line[0] == "#"):
                continue
            vcf_cut.write(line + "\n")
        vcf_cut.close()
        vcf_df = pd.read_csv("{0}_abbrev.txt".format(filename), sep=',', header=0)
        return vcf_df
    
    # summary: get a pileup read counts dataframe
    # parameter: read counts file
    # returns: a DataFrame with all the read counts information
    def getReadCounts(self):
        readCounts = self.readCounts
        df = pd.read_csv(readCounts, sep=",")
        return df

    # summary: get complimentary base
    # parameter: character base pair
    # returns: complimentary base pair
    def comp(self, char):
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
    def ref_alt(self, df, vcf_df, peakLocs):
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
            key = vcf_df.iloc[[i]]["POS"].item()
            snpIndex = peakLocs[vcf_df.iloc[[i]]["CHROM"].item()][key]
            # if this location does not fall in a peak none of this matters
            if snpIndex == 0:
                continue
            # set refseq and altseq to the items in each if they are an alphanumeric characte
            refseq = "".join(c for c in vcf_df.iloc[[i]]["REF"].item() if c.isalnum())
            altseq = "".join(c for c in vcf_df.iloc[[i]]["ALT"].item() if c.isalnum())
            # indel events actually have strings of longer than 1, therefore just continue and keep track
            if len(refseq) > 1 or len(altseq) > 1:
                indelCount += 1
                continue
            bases = set(['a', 'g', 'c', 't'])
            if (refseq.lower() not in bases) or (altseq.lower() not in bases):
                continue
            # reference and alternate counts are gathered by using this dictionary referencing method
            ref = df.iloc[i][refseq.upper()].item() + df.iloc[i][self.comp(refseq).lower()].item()
            alt = df.iloc[i][altseq.upper()].item() + df.iloc[i][self.comp(altseq).lower()].item()
            # if the counts are both 0, then clearly there is an issue (snp with no read counts)
            if ref == 0 and alt == 0:
                err.append([ref, alt])
                continue
            # if sum of counts is less than 3 then skip due to lower power
            if ref + alt < 3:
                continue
            # if the vcf_index is not a key then create a list
            if snpIndex not in ref_alt.keys():
                ref_alt[snpIndex] = []
            # append the read counts tuple to the list
            ref_alt[snpIndex].append((ref, alt, key))
        return ref_alt, err

    # summary: generates a DataFrame with the start and end of the peaks from the peaks file
    # parameters: file of the peaks formatted: "CHR START END ID SAMPLE1 ..." tab separated, and the delimiter
    # returns: generated DataFrame
    def gen_peaks_df(self, peaksFile, sep=","):
        peaks = pd.read_csv(peaksFile, sep=sep)
        size = peaks.shape[0]
        d = {}
        for i in range(size):
            d[i] = [peaks.iloc[[i]]["CHR"].item(), str(int(peaks.iloc[[i]]["START"].item()) - 150), str(int(peaks.iloc[[i]]["END"].item()) + 150)]
        peaks_df = pd.DataFrame(list(d.values()), columns=["CHR", "START", "END"])
        return peaks_df

    # summary: makes a dictionary of the snps that fall under a peak (numpy zeros with values only stored if there is a snp there)
    # parameters: the peaks DataFrame
    # returns: dictionary of the format {chr: np.array()}
    def peaksArray(self, peaks_df):
        # dictionary to store all the outputs
        d = {}
        # index is the index of snps
        index = 0
        # iterate through the peaks and see if any of the snps fall into them
        peak_num = 0
        last = peaks_df.shape[0]
        mx = int(max([int(i) for i in peaks_df["END"].values]))
        for i in peaks_df["CHR"].tolist():
            d[i] = np.zeros(mx+1)
        for i in range(last):
            start = int(peaks_df.iloc[[i]]["START"].item())
            end = int(peaks_df.iloc[[i]]["END"].item())
            d[peaks_df.iloc[[i]]["CHR"][i]][start:end] = i
        return d

    # summary: get the counts for each snp that falls under peaks
    # paramters: reference/alternate dictionary, location of peaks {chr: np.array()}
    # returns: dictionary of the formate {chr: {peakIndex: [(ref, alt), ...], ...}, ...}
    def peak_counts(ref_alt_dict, peak_locs):
        refseq = ""
        altseq = ""
        alt = 0
        ref = 0
        ase = 0
        # finalDict - final dictionary that is returned
        finalDict = {}
        # n is the number of snps with a valid number of ref and alt counts (not 0 for both)
        n = len(ref_alt_dict.keys())
        # ignore indels
        for chrom in ref_alt_dict.keys():
            # peaks_counts format {index of peak: [(ref, alt), ...]}
            peak_counts = {}
            for index in ref_alt_dict[chrom].keys():
                # make a list if it doesn't exist
                if peak_locs[chrom][index] == 0:
                    continue
                # if it does exist
                else:
                    # if the peak number is not 
                    if int(peak_locs[chrom][index]) not in peak_counts.keys():
                        peak_counts[int(peak_locs[chrom][index])] = []
                    peak_counts[peak_locs[chrom][index]].extend(ref_alt_dict[chrom][index])
            finalDict[chrom] = peak_counts
        return finalDict
       
    # summary: get the observed ase at each snp site in the reference/alternate dictionary
    # parameters: reference/alternate dictonary
    # returns: dictionary with the observed ASE values
    def getObservedASE(self, ref_alt):
        ase_o = {}
        ase = 0
        for i in ref_alt.keys():
            ase_o[i] = []
            for counts in ref_alt[i]:
                if np.NaN in counts:
                    continue
                ref, alt, pos = counts
                ase = abs(ref - alt) / float(ref + alt)
                ase_o[i].append(ase)
        return ase_o
        
    # summary: get the binomial simulation of ase given the peak counts
    # parameters: reference/alternate dictionary {peakIndex: [refCount, altCount, POS]}, probability (default: 1/2), trials (default: 10000)
    # returns: binomial ase simulated
    def getSimulatedASE(self, ref_alt, p=1/2, trials=10000):
        all_counts = 0
        # xbr => reference binomial simulated
        # xba => alternate binomial simulated
        xbr = 0
        xba = 0
        # ase_b is binomial ase
        ase_b = {}
        for key in ref_alt.keys():
            ase_b[key] = []
            ase_b_arr = np.zeros(trials)
            num_snps = float(len(ref_alt[key]))
            for snp in ref_alt[key]:
                if np.NaN in snp:
                    continue
                all_counts = snp[0] + snp[1]
                xbr = np.random.binomial(all_counts, p, trials)
                xba = np.ones(trials) * float(all_counts) - xbr
                ase_b_arr += abs(xbr - xba)/float(all_counts)
            ase_b_arr = ase_b_arr / num_snps
            ase_b[key] = ase_b_arr
        return ase_b
        
    # summary: get the p values of ase
    def getPVals(self, trials, p, ref_alt, ase_b_peak, ase_o_peak, chrom):
        p_vals = []
        for peakIndex in ref_alt:
            p_value = 0.0
            if peakIndex in ase_b_peak and peakIndex in ase_o_peak:
                p_value = float(sum([np.sum(ase_b_peak[peakIndex][i] >= ase_o_peak[peakIndex][0]) for i in range(len(ase_b_peak[peakIndex]))])) + 1.0
            p_value /= float(trials + 1)
            p_vals.append((peakIndex, p_value))
        return p_vals

    # runs the main method
    def run(self):
        # get DataFrame of the read counts piled up at each het site
        # columns of df
        # chr loc ref A	T C G a t c g Insertion Deletion
        getReadCountsStart = timeit.default_timer()
        print("Get Read Counts")
        pileupDf = self.getReadCounts()
        getReadCountsEnd = timeit.default_timer()
        print("pileupDF shape: {0}".format(str(pileupDf.shape)))
        print("time for Get Read Counts: {0}".format(str(getReadCountsEnd - getReadCountsStart)))

        # get a dataframe containing the vcf information
        # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT MSSM_13 ...
        processVCFStart = timeit.default_timer()
        print("Process VCF")
        vcf_df = self.proc_vcf(self.vcf)
        processVCFEnd = timeit.default_timer()
        print("vcf_df shape: {0}".format(str(vcf_df.shape)))
        print("time for Process VCF: {0}".format(str(processVCFEnd - processVCFStart)))
        
        # get peaks from the peaks file
        # CHR START END
        genPeaksDFStart = timeit.default_timer()
        print("Generate Peaks DF")
        peaks = self.gen_peaks_df(self.peaks, ",")
        genPeaksDFEnd = timeit.default_timer()
        print("peaks shape: {0}".format(str(peaks.shape)))
        print("time for Generate Peaks DF: {0}".format(str(genPeaksDFEnd - genPeaksDFStart)))

        # get dictionary storing the numpy array of where all the peaks are
        # {CHR: np.array(location of the snps)}
        peakLocsStart = timeit.default_timer()
        print("Get Dictionary Storing Numpy Array of Peaks")
        peakLocs = self.peaksArray(peaks)
        peakLocsEnd = timeit.default_timer()
        print("peakLocs len: {0}".format(str(len(peakLocs))))
        print("time for Get Dictionary Storing Numpy Array of Peaks: {0}".format(str(peakLocsEnd - peakLocsStart)))

        # get a dictionary of the format
        # filter for all the peak locations within this method
        # {CHR: {peakIndex: [(ref, alt, pos), ...]},  ...}
        filterPeakLocsStart = timeit.default_timer()
        print("Filter All The Peak Locations")
        ref_alt_dict = {}
        err_dict = {}
        for chrom in peakLocs.keys():
            input_df = pileupDf.loc[pileupDf['chr'] == chrom]
            input_vcf_df = vcf_df.loc[vcf_df['CHROM'] == chrom]
            ref_alt, err = self.ref_alt(input_df, input_vcf_df, peakLocs)
            ref_alt_dict[chrom] = ref_alt
            err_dict[chrom] = err
        # pd.DataFrame.from_dict(data=ref_alt_dict, orient='index').to_csv("{0}/ref_alt_dict.csv".format(self.outDir), header=False)
        ref_alt_df = pd.DataFrame.from_dict(data=ref_alt_dict, orient='index')
        # remove all NaN only columns
        ref_alt_df = ref_alt_df.dropna(axis=1, how='all')
        # transpose to make the columns the rows and vice versa
        ref_alt_df = ref_alt_df.T
        ref_alt_df.to_json("{0}/ref_alt_dict.json".format(self.outDir), orient="split")
        filterPeakLocsEnd = timeit.default_timer()
        print("ref_alt_dict len: {0}".format(str(len(ref_alt_dict))))
        print("time for Filter All The Peak Locations: {0}".format(str(filterPeakLocsEnd - filterPeakLocsStart)))

        # get the snps under a peak
        # {CHR: {peakIndex: [(ref, alt), ...]} ...}
        # peakCounts = self.peak_counts(ref_alt_dict, peakLocs)

        # get the observed ase from the samples on a snp-by-snp basis
        # {CHR: {peakIndex: [ase_o, ...], ...} ...}
        getObservedASEStart = timeit.default_timer()
        print("Get Observed ASE")
        ase_o_dict = {}
        for chrom in peakLocs.keys():
            ase_o = self.getObservedASE(ref_alt_dict[chrom])
            ase_o_dict[chrom] = ase_o
        getObservedASEEnd= timeit.default_timer()
        pd.DataFrame.from_dict(data=ase_o_dict, orient='index').to_json("{0}/ase_o_dict{1}.json".format(self.outDir,self.chrom), orient="split")
        print("ase_o_dict len: {0}".format(str(len(ase_o_dict))))
        print("time for Get Observed ASE: {0}".format(str(getObservedASEEnd - getObservedASEStart)))
        
        # get the ase overall per peak
        # {CHR: {peakIndex: ase_o_overall, ...} ...}
        # ase_o_peak_dict = {}
        # for chrom in ase_o_dict.keys():
        #    ase_o_peak = {}
        #    for peakIndex in ase_o_dict[chrom].keys():
        #        ase_o_peak[peakIndex] = sum(ase_o_dict[chrom][peakIndex])/len(ase_o_dict[chrom][peakIndex])
        #    ase_o_peak_dict[chrom] = ase_o_peak

        
        # WORK ON GETTING ASE_B WORKING IN THEORY
        # get the simulated ase overall per peak
        getSimulatedASEStart = timeit.default_timer()
        trials = 100000
        p = 0.5
        print("Get Simulated ASE")
        ase_b_dict = {}
        # {CHR: {peakIndex: ase_b_overall, ...} ...}
        for chrom in peakLocs.keys():
            ase_b = self.getSimulatedASE(ref_alt_dict[chrom], p, trials)
            ase_b_dict[chrom] = ase_b
        getSimulatedASEEnd = timeit.default_timer()
        pd.DataFrame.from_dict(data=ase_b_dict, orient='index').to_json("{0}/ase_b_dict{1}.json".format(self.outDir,self.chrom), orient="split")
        print("ase_b_dict len: {0}".format(str(len(ase_b_dict))))
        print("time for Get Simulated ASE: {0}".format(str(getSimulatedASEEnd - getSimulatedASEStart)))

        # get p values
        getPValuesStart = timeit.default_timer()
        print("Get P Values")
        pVals = {}
        for chrom in peakLocs.keys():
            pVals[chrom] = self.getPVals(trials, p, ref_alt_dict[chrom], ase_b_dict[chrom], ase_o_dict[chrom], chrom)
        pd.DataFrame.from_dict(data=pVals, orient='index').to_json("{0}/p_values{1}.json".format(self.outDir,self.chrom), orient="split")
        getPValuesEnd = timeit.default_timer()
        print("pVals shape: {0}".format(len(pVals)))
        print("time for Get P Values: {0}".format(str(getPValuesEnd - getPValuesStart)))
        print("Done")

aseModel = ase_simpler()
aseModel.run()
