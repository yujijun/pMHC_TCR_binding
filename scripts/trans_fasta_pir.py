#!/usr/bin/python3
from modeller import *
import os
env = environ()
# 将一个FASTA序列转化为PIR格式的序列;
# 如果文件中有多条序列，输出结果将是比对后的结果;
# 一般用于一条序列的格式转化
input_fasta = input("Input FASTA sequence file:")
input_name = input_fasta.split(".")[-2]
aln = alignment(env)
aln = alignment(env, file=input_fasta, alignment_format='FASTA')
aln.write(file=input_name+'.ali', alignment_format='PIR')

# os.system("sed '1d' input_name+'.ali'")
