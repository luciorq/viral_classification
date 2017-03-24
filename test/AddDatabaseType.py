import sys
import os
import csv

file_list = os.listdir("data/")
Type = {}

Virus_fasta = "data/Drosophila_viruses.fasta"
with open(Virus_fasta, "r") as fp:
	content = fp.read()
for file in file_list:
	if not file.endswith('.fsa'):
		continue
	seq_id = file.split('.')
	if seq_id[0] in content:
		Type[seq_id[0]] = "Virus"

EVE_fasta = "data/EVEs.fasta"
with open(EVE_fasta, "r") as fp:
	content = fp.read()
for file in file_list:
	if not file.endswith('.fsa'):
		continue
	seq_id = file.split('.')
	if seq_id[0] in content:
		Type[seq_id[0]] = "EVE"

TE_fasta = "data/dmel-all-transposon-r6.13.fasta"
with open(TE_fasta, "r") as fp:
	content = fp.read()
for file in file_list:
	if not file.endswith('.fsa'):
		continue
	seq_id = file.split('.')
	if seq_id[0] in content:
		Type[seq_id[0]] = "TE"

TE2_fasta = "data/Aaegypti.TEs_TEfam.fasta"
with open(TE2_fasta, "r") as fp:
	content = fp.read()
for file in file_list:
	if not file.endswith('.fsa'):
		continue
	seq_id = file.split('.')
	if seq_id[0] in content:
		Type[seq_id[0]] = "TE"

database = "data/db1.db"
with open(database, "r") as fp_db:

	db = csv.reader(fp_db, delimiter= "\t")

	db = zip(*db)
	with open("data/db1_decision_tree.db", "w") as fp_tree:
		for row in db:
			actual_row = list(row)
			appended = 0
			if actual_row[0] == "lib_name":
				actual_row.append("Type")
				appended = 1
			for name in Type:
				if name in actual_row[0]:
					#print Type[name]
					actual_row.append(Type[name])
					appended = 1	
			if appended == 0:
				actual_row.append("Unknown")
			for item in actual_row:
				fp_tree.write("%s\t" % item)
			fp_tree.write("\n")
		


