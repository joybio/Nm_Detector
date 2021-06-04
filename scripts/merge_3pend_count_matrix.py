#!/root/miniconda3/bin/python

__date__ = "2021-5-20"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import sys
import os
from sys import argv
import collections
from collections import defaultdict

result_path = argv[1].split("/")
print(result_path)
result_dir = '/'.join(result_path[:-1])
print(result_dir)
result_file = result_dir + "/3pend_count_matrix.xls"
result_file2 = result_dir + "/3pend_counts_matrix.xls"
out = open(result_file,"w")
out2 = open(result_file2,"w")

#write head in output files
out.write("chr\tstart\tstop")
for i in range(1,len(argv)):
	sample = os.path.basename(argv[i]).split(".")
	out.write("\t" + sample[0])
out.write("\n")

out2.write("position")
for i in range(1,len(argv)):
	sample = os.path.basename(argv[i]).split(".")
	out2.write("\t" + sample[0])
out2.write("\n")

#extract all 3pend sites from all samples
pos_set = set()
for i in range(1,len(argv)):
	data = open(argv[i],'r')
	for i in data:
		i = i.strip().split("\t")
		key = i[0]+"\t" +i[1]+"\t" +i[2]
		pos_set.add(key)
	data.close()
print("Total reads end number:%d"%len(pos_set))
data_dict = defaultdict(list)
for i in range(1,len(argv)):
	data = open(argv[i],'r')
	if i == 1:
		for j in data:
			j = j.strip().split("\t")
			key = str(j[0]+"\t" +j[1]+"\t" +j[2])
			value = str(j[3])
			data_dict[key].append(value)
		for k in pos_set:
			if k not in data_dict.keys():
				data_dict[k].append("0")
			else:
				pass
		
	else:
		sub_data_dict = {}
		for j in data:
			j = j.strip().split("\t")
			key = str(j[0]+"\t" +j[1]+"\t" +j[2])
			value = str(j[3])
			sub_data_dict[key] = value
			data_dict[key].append(value)
		for k in pos_set:
			if k not in sub_data_dict.keys():
				data_dict[k].append("0")
			
	data.close()
		
for i in data_dict.keys():
	reads_end_cnt = "\t".join(data_dict[i])
	pos = i.replace("\t","__")
	out.write(i + "\t" + reads_end_cnt + "\n")
	out2.write(pos + "\t" + reads_end_cnt + "\n")
out.close()
out2.close()
sys.exit()
