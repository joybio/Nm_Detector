#!/root/miniconda3/bin/python

"""
extract known stop sites: stop_codon; transcript_stop_sites; start_codon; tRNA_stop_sites; other_RNA_stop_sites.
start_codon; stop_codon; transcript_stop_sites, other_RNA_stop_sites: gtf file
tRNA_start_sites; tRNA_stop_site: http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa
"""

__date__ = "2021-5-24"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import re
import os
import optparse
from optparse import OptionParser
parser = OptionParser('Usage: %prog -b input -g gtf_file -t tRNA_file -o output')
parser.add_option('-b','--bed',
		dest='bed',
		help='pre_results in bed format')

parser.add_option('-g','--gtf',
		dest='gtf',
		help='gtf file')

parser.add_option('-t','--tRNA',
		dest='tRNA',
		help='tRNA file:hg38-mature-tRNAs.fa')

parser.add_option('-o','--out',
		dest='out',
		help='out annotation file')
(options,args) = parser.parse_args()

#extract known stop sites.
#start_codon; stop_codon; transcript_start_sites; transcript_stop_sites
gtf = open(options.gtf,"r")
known_sites = set()
for i in gtf:
	if i.startswith("#"):
		pass
	else:
		i=i.strip().split("\t")
		i[3] = int(i[3])
		i[4] = int(i[4])
		if i[2] == "gene" or i[2] == "transcript" or i[2] == "start_codon" or \
		i[2] == "stop_codon" or i[2] == "five_prime_utr" or i[2] == "three_prime_utr":
			if i[6] == "+":
				#start
				known_site_1 = i[0] + "\t" + str(i[3]-1) + "\t" + str(i[3])
				known_site_2 = i[0] + "\t" + str(i[3]-2) + "\t" + str(i[3]-1)
				known_site_3 = i[0] + "\t" + str(i[3]) + "\t" + str(i[3]+1)
				known_site_11 = i[0] + "\t" + str(i[3]+1) + "\t" + str(i[3]+2)
				known_site_12 = i[0] + "\t" + str(i[3]+2) + "\t" + str(i[3]+3)
				known_site_13 = i[0] + "\t" + str(i[3]-1) + "\t" + str(i[3]-2)
				known_site_14 = i[0] + "\t" + str(i[3]-2) + "\t" + str(i[3]-3)
				#stop
				known_site_4 = i[0] + "\t" + str(i[4]-1) + "\t" + str(i[4])
				known_site_5 = i[0] + "\t" + str(i[4]-2) + "\t" + str(i[4]-1)
				known_site_6 = i[0] + "\t" + str(i[4]) + "\t" + str(i[4]+1)
				known_site_7 = i[0] + "\t" + str(i[4]-3) + "\t" + str(i[4]-2)
				known_site_8 = i[0] + "\t" + str(i[4]-4) + "\t" + str(i[4]-3)
				known_site_9 = i[0] + "\t" + str(i[4]+1) + "\t" + str(i[4]+2)
				known_site_10 = i[0] + "\t" + str(i[4]+2) + "\t" + str(i[4]+3)
				known_sites.add(known_site_1)
				known_sites.add(known_site_2)
				known_sites.add(known_site_3)
				known_sites.add(known_site_4)
				known_sites.add(known_site_5)
				known_sites.add(known_site_6)
				known_sites.add(known_site_7)
				known_sites.add(known_site_8)
				known_sites.add(known_site_9)
				known_sites.add(known_site_10)
				known_sites.add(known_site_11)
				known_sites.add(known_site_12)
				known_sites.add(known_site_13)
				known_sites.add(known_site_14)
gtf.close()
tRNA = open(options.tRNA,"r")
for i in tRNA:
	if i.startswith(">"):
		i = i.strip().split("chr")
		j = i[2].split(" (")
#		print(j)
		m = j[0].split(":")
		chrom = m[0]
		n = m[1].split("-")
		start = int(n[0])
		stop = int(n[1])
		known_site_1 = chrom + "\t" + str(start-1) + "\t" + str(start)
		known_site_2 = chrom + "\t" + str(start-2) + "\t" + str(start-1) 
		known_site_3 = chrom + "\t" + str(start) + "\t" + str(start+1)
		known_site_11 = chrom + "\t" + str(start+1) + "\t" + str(start+2)
		known_site_12 = chrom + "\t" + str(start+2) + "\t" + str(start+3)
		known_site_13 = chrom + "\t" + str(start-1) + "\t" + str(start-2)
		known_site_14 = chrom + "\t" + str(start-2) + "\t" + str(start-3)
		#stop
		known_site_4 = chrom + "\t" + str(stop-1) + "\t" + str(stop)
		known_site_5 = chrom + "\t" + str(stop-2) + "\t" + str(stop-1)
		known_site_6 = chrom + "\t" + str(stop) + "\t" + str(stop+1)
		known_site_7 = chrom + "\t" + str(stop+1) + "\t" + str(stop+2)
		known_site_8 = chrom + "\t" + str(stop+2) + "\t" + str(stop+3)
		known_site_9 = chrom + "\t" + str(stop-3) + "\t" + str(stop-2)
		known_site_10 = chrom + "\t" + str(stop-4) + "\t" + str(stop-3)
		known_sites.add(known_site_1)
		known_sites.add(known_site_2)
		known_sites.add(known_site_3)
		known_sites.add(known_site_4)
		known_sites.add(known_site_5)
		known_sites.add(known_site_6)
		known_sites.add(known_site_7)
		known_sites.add(known_site_8)
		known_sites.add(known_site_9)
		known_sites.add(known_site_10)
		known_sites.add(known_site_11)
		known_sites.add(known_site_12)
		known_sites.add(known_site_13)
		known_sites.add(known_site_14)
tRNA.close()

data = open(options.bed,"r")
out = open(options.out,"w")
head = data.readline()
out.write(head)
for i in data:
	line = i
	i = i.strip().split("\t")
	pos = i[10] + "\t" + str(i[11]) + "\t" + str(i[12])
	if pos in known_sites:
		pass
	else:
		out.write(line)

data.close()
out.close()

