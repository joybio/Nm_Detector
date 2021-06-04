#!/root/miniconda3/bin/python

__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#import numpy as np
#import pandas as pd
import collections
from collections import defaultdict
import optparse
from optparse import OptionParser

parser = OptionParser("Usage: %prog ")
parser.add_option("-i","--input",dest = "input",
		help = "Input file; the mapping results of reads1. If hisat2 is used, use param '-k 1' and grep NH:i:1 in the sam file.")
parser.add_option("-o","--output",dest = "out",
		help = "Output file (bed); the results of 3pend site.")
parser.add_option("-f","--file",dest = "file",
		help = "Output file (format bed); the results of 3pend site.")
(options,args) = parser.parse_args()

data = open(options.input,"r")
out = open(options.out,"w")
formatfile = open(options.file,"w")
#extract 3pend based on the strandness.
pos_dict = defaultdict(int)

for i in data:
	i = i.strip().split("\t")
	position = i[0]+ "\t" + i[1] +"\t"+ i[2]
	pos_dict[position] += 1
for i in pos_dict.keys():
	out.write(i + "\t" + str(pos_dict[i]) + "\n")
for i in pos_dict.keys():
	formatfile.write("chr"+i + "\t" + str(pos_dict[i]) + "\n")
data.close()
out.close()
formatfile.close()












