###################################
#    统计Allelelist.txt中出现频率	  #
#	    排在前十的等位基因			  #
###################################
fr = open("Allelelist-2021-01-18.txt","r")

allele_DPA1_statistic = dict()
allele_DPB1_statistic = dict()
allele_DQA1_statistic = dict()
allele_DQA2_statistic = dict()
allele_DQB1_statistic = dict()
allele_DRA_statistic = dict()
allele_DRB1_statistic = dict()

for line in fr.readlines():
	if line.startswith('#'):
		continue
	elif line.startswith('AlleleID'):
		continue

	line = line.strip()
	line_split = line.split(',')
	allele_loc = line_split[1].split('*')[0]
	allele_name_split = line_split[1].split(':')
	allele_name_2 = ':'.join(allele_name_split[0:2])

	if allele_loc == 'DPA1':
		if allele_name_2 in allele_DPA1_statistic.keys():
			allele_DPA1_statistic[allele_name_2] += 1
		else:
			allele_DPA1_statistic[allele_name_2] = 1
	elif allele_loc == 'DPB1':
		if allele_name_2 in allele_DPB1_statistic.keys():
			allele_DPB1_statistic[allele_name_2] += 1
		else:
			allele_DPB1_statistic[allele_name_2] = 1
	elif allele_loc == 'DQA1':
		if allele_name_2 in allele_DQA1_statistic.keys():
			allele_DQA1_statistic[allele_name_2] += 1
		else:
			allele_DQA1_statistic[allele_name_2] = 1
	elif allele_loc == 'DQA2':
		if allele_name_2 in allele_DQA2_statistic.keys():
			allele_DQA2_statistic[allele_name_2] += 1
		else:
			allele_DQA2_statistic[allele_name_2] = 1
	elif allele_loc == 'DQB1':
		if allele_name_2 in allele_DQB1_statistic.keys():
			allele_DQB1_statistic[allele_name_2] += 1
		else:
			allele_DQB1_statistic[allele_name_2] = 1
	elif allele_loc == 'DRA':
		if allele_name_2 in allele_DRA_statistic.keys():
			allele_DRA_statistic[allele_name_2] += 1
		else:
			allele_DRA_statistic[allele_name_2] = 1
	elif allele_loc == 'DRB1':
		if allele_name_2 in allele_DRB1_statistic.keys():
			allele_DRB1_statistic[allele_name_2] += 1
		else:
			allele_DRB1_statistic[allele_name_2] = 1
fr.close()

def my_sort(allele_statistic):
	allele_statistic = sorted(allele_statistic.items(), key=lambda x:x[1], reverse=True)
	allele_statistic = dict(allele_statistic[:10])
	allele_statistic = sorted(allele_statistic.items(), key=lambda x:x[0], reverse=False)
	return allele_statistic

allele_DPA1_statistic = my_sort(allele_DPA1_statistic)
allele_DPB1_statistic = my_sort(allele_DPB1_statistic)
allele_DQA1_statistic = my_sort(allele_DQA1_statistic)
allele_DQA2_statistic = my_sort(allele_DQA2_statistic)
allele_DQB1_statistic = my_sort(allele_DQB1_statistic)
allele_DRA_statistic = my_sort(allele_DRA_statistic)
allele_DRB1_statistic = my_sort(allele_DRB1_statistic)

allele_statistic = allele_DPA1_statistic + allele_DPB1_statistic \
					+ allele_DQA1_statistic + allele_DQA2_statistic + allele_DQB1_statistic \
					+ allele_DRA_statistic + allele_DRB1_statistic \

allele_seq = dict()
flag = 0
for i in allele_statistic:
	allele_name = i[0]
	allele_seq[allele_name] = ''
	fr = open("hla_prot.fasta","r")
	for line in fr.readlines():
		line = line.strip()
		if line.startswith('>'):
			line_split = line.split(' ')
			allele_name_split = line_split[1].split(":")
			seq_length = line.split(' ')[2]
			if ':'.join(allele_name_split[0:2]) == allele_name:
				if int(seq_length) > len(allele_seq[allele_name]):
					allele_seq[allele_name] = ">"+allele_name + " " + seq_length+" bp,"
					flag = 1
				else:
					flag = 0
			else:
				flag = 0

		elif(flag == 1):
			allele_seq[allele_name] += line

	fr.close()

fw = open("Top_MHCII_alleles.csv","w")
for key in allele_seq.keys():
	fw.write("HLA-"+key+","+allele_seq[key]+"\n")

fw.close()

