# this file removes all the blacklist regions
# a peak will be removed as long as it has intersection with the blacklist area.
# Input: results/peak/peakFile, data/blackList_mm9/blacklists_chr6

import os

peakFile = "../results/peaks/peaks_00h/00h_peaks.bed"

#os.chdir("../../Documents/BB2289_projectInMolecularBio/scripts")

# Read blacklist file
with open("../data/blackList_mm9/mm9_chr6_blacklist.bed",'r') as f:
    blackList = f.read().splitlines() 
    blackRegions = []
    for lines in blackList:
        blackRegions.append(lines.split('\t'))
    f.close()

# Read peakfile
with open(peakFile,'r') as g:
    peakList = g.read().splitlines() 
    peaks = []
    for lines in peakList:
        peaks.append(lines.split('\t'))
    f.close()

# Find black regions in peaks
newPeaks = []
for peak in peaks:
    for black in blackRegions:
        if (int(peak[1]) > int(black[1]) and int(peak[1]) < int(black[2])):
            pass
        if (int(peak[2]) > int(black[1]) and int(peak[2]) < int(black[2])):
            pass
    
    #print(peak[2])
        
    