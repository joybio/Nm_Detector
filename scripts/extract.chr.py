#!/root/miniconda3/bin/python

"""
chromsome size: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
bedGraphToBigwig: https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/bedGraphToBigWig
mapping reference is different with chromsome size, before use bedGraphToBigwig, you should format your bed file to chromsome size. 
"""

__date__ = "2021-5-24"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-b','--bed',
                dest='bed',
                help='bed file. Here means the reads end coverage.')

parser.add_option('-o','--out',
                dest='out',
                help='out file.')

(options,args) = parser.parse_args()

data = open("/home/l/backup1/joybio/script/bed2bw/hg38.chrom.sizes","r")
size_set = set()
for i in data:
	i = i.strip().split("\t")
	size_set.add(i[0])
data.close()
bedfile = open(options.bed,"r")
outfile = open(options.out,"w")
bed_dict = {}
for i in bedfile:
	line = i
	i = i.strip().split("\t")
	key = i[0]
	if key in size_set:
		outfile.write(line)
bedfile.close()
outfile.close()



