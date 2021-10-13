#!/usr/bin/python3
from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A'],
                             renumber_residues=[1])
        # Another way to label individual chains:
        # self.chains[0].name = 'A'
        # self.chains[1].name = 'B'

env = environ()


aln = alignment(env)
code =[]
aln_file = input("alignment sequence file(.ali):")
read_aln = modfile.File(aln_file, 'r')
while aln.read_one(read_aln, alignment_format='PIR'):
	code.append(aln[0].code)
	
read_aln.close()
template_code = code[0]
target_code = code[1]

a = MyModel(env, 
			  alnfile=aln_file,
              knowns=template_code, 
              sequence=target_code,
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
