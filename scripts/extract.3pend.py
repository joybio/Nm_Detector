#!/root/miniconda3/bin/python

"""
1. The orientation of reads1 should be same with the gene orientation, thus, we can use the 3pend of the reads1 to stand for the stop site, which should be the Nm site.
2. the length of the reads1 is very short (normaly 20-40), thereby, reads1 and reads2 stand for the same fragement. No need to calculate the 3pend of reads2. 
"""

__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser

parser = OptionParser("Usage: %prog ")
parser.add_option("-i","--input",dest = "input",
		help = "Input file; the mapping results of reads1. If hisat2 is used, use param '-k 1' and grep NH:i:1 in the sam file.")
parser.add_option("-o","--output",dest = "out",
		help = "Output file; the results of 3pend site.")

(options,args) = parser.parse_args()

data = open(options.input,"r")
out = open(options.out,"w")

#extract 3pend based on the strandness.

for i in data:
	i = i.strip().split("\t")
	reads_seq = i[9]
	chromosme = i[2]
	start = int(i[3])
	stop = int(i[3]) + len(reads_seq)
	# strandness == "+";
	#for example: sam file: start=2414 seq=NTGTTGTGAGTTGCAGAGACAGGAAAG。
	#bedtools getfasta results: TGTTGTGAGTTGCAGAGACAGGAAAGATAGAAGACGTTCCAT
	#start nuclertide 2415. pos=2415-2416
	#for example: extract region 1-15, results: TGTCGTCCGTTTCA，not:ATGTCGTCCGTTTCA
	#for example: bedtools extract region 1-28 results: TGTCGTCCGTTTCAGGAAATGCCTCTT, total len is 27bp。
	if "XS:A:+" in i:
		out.write(chromosme + "\t" + str(stop) + "\t" + str(stop+1) + "\n")
	elif "XS:A:-" in i:
		out.write(chromosme + "\t" + str(start-1) + "\t" + str(start) + "\n")
data.close()
out.close()













