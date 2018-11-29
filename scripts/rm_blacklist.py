# this file removes all the blacklist regions
# a peak will be removed as long as it has intersection with the blacklist area.
# Input: results/peak/peakFile, data/blackList_mm9/blacklists_chr6

import os
os.path.expanduser("~/Documents/BB2289_projectInMolecularBio")
os.getcwd()

peakFile = "../results/peaks/peaks_12h/12h_peaks.bed"
outFile = "../results/peaks/peaks_12h/12h_peaks_true.bed"

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

# Instead of finding black regions, create a white list
first = int(peaks[0][2])-1
last = int(peaks[len(peaks)-1][2])+1

whiteStart= []
whiteStart.append(first)
whiteEnd = []
for region in blackRegions:
    whiteStart.append(int(region[2]))
    whiteEnd.append(int(region[1]))
whiteEnd.append(last)

# Find true peaks
count = 0
peakTrue = []
for peak in peaks:
    for i in range(0,len(whiteEnd)):
        if int(peak[1])>whiteStart[i] and int(peak[2])<whiteEnd[i]:
            count = count + 1
            peakTrue.append(peak)
            
print("the number of valid peaks: "+str(count))

# Save as bed file
with open(outFile, 'w') as f:
    for peak in peakTrue:
        for string in peak:
            f.write(str(string) + '\t')
        f.write('\n')
f.close()