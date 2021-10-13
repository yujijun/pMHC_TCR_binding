###################################
#    统计Allelelist.txt中出现频率	  #
#	    排在前十的等位基因			  #
###################################
fr = open("Allelelist-2021-01-18.txt","r")

allele_A_statistic = dict()
allele_B_statistic = dict()
allele_C_statistic = dict()

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

	if allele_loc == 'A':
		if allele_name_2 in allele_A_statistic.keys():
			allele_A_statistic[allele_name_2] += 1
		else:
			allele_A_statistic[allele_name_2] = 1
	elif allele_loc == 'B':
		if allele_name_2 in allele_B_statistic.keys():
			allele_B_statistic[allele_name_2] += 1
		else:
			allele_B_statistic[allele_name_2] = 1
	elif allele_loc == 'C':
		if allele_name_2 in allele_C_statistic.keys():
			allele_C_statistic[allele_name_2] += 1
		else:
			allele_C_statistic[allele_name_2] = 1
fr.close()

def my_sort(allele_statistic):
	allele_statistic = sorted(allele_statistic.items(), key=lambda x:x[1], reverse=True)
	allele_statistic = dict(allele_statistic[:10])
	allele_statistic = sorted(allele_statistic.items(), key=lambda x:x[0], reverse=False)
	return allele_statistic

allele_A_statistic = my_sort(allele_A_statistic)
allele_B_statistic = my_sort(allele_B_statistic)
allele_C_statistic = my_sort(allele_C_statistic)

allele_statistic = allele_A_statistic + allele_B_statistic + allele_C_statistic

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

fw = open("Top_MHCI_alleles.csv","w")
for key in allele_seq.keys():
	fw.write("HLA-"+key+","+allele_seq[key]+"\n")

fw.close()

