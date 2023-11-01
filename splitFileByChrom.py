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

class split():
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-rc", "--readCounts", dest="readCounts", help="get read counts from the read counts file")
        parser.add_argument("-v", "--vcf", dest="vcf", help="get het site information using the vcf")
        parser.add_argument("-p", "--peaks", dest="peaks", help="get peaks with a start and end index")
        parser.add_argument("-o", "--outputDirectory", dest="outDir", help="output directory for pVals dictionary (outputs csv)")
        args = parser.parse_args()
        self.readCounts = args.readCounts
        self.vcf = args.vcf
        self.peaks = args.peaks
        self.outDir = args.outDir

    def splitReadCounts(self):
        df = pd.read_csv(self.readCounts)
        chrms = df['chr'].unique()
        for chrm in chrms:
            outputDf = df.loc[df['chr'] == chrm]
            outputDf.to_csv(self.outDir + "/output/{0}.txt".format(chrm))
    
    def splitVCF(self):
        df = pd.read_csv(self.vcf)
        chrms = df['CHROM'].unique()
        for chrm in chrms:
            outputDf = df.loc[df['CHROM'] == chrm]
            outputDf.to_csv(self.outDir + "/vcf/{0}.txt".format(chrm))
        
    def splitPeaks(self):
        df = pd.read_csv(self.peaks)
        chrms = df['CHR'].unique()
        for chrm in chrms:
            outputDf = df.loc[df['CHR'] == chrm]
            outputDf.to_csv(self.outDir + "/peaks/{0}.txt".format(chrm))

    def run(self):
        self.splitReadCounts()
        self.splitVCF()
        self.splitPeaks()

obj = split()
obj.run()
