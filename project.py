#!/usr/bin/env python3
from __future__ import division

import gffutils
import os, itertools, re

def isheader(line):
	return line[0] == '>'

def aspairs(f):
	seq_id = ''
	sequence = ''
	for header,group in itertools.groupby(f, isheader):
		if header:
			line = next(group)
			seq_id = line[1:].split()[0]
		else:
			sequence = ''.join(line.strip() for line in group)
			yield seq_id, sequence

#Make a database from the chr1_gff3.gff3 file:
#force-True will rewrite the db
#merge_strategy="merge" will merge duplicate data
#from_string=True shows that the data is the actual data to use

data = ''
with open('chr1_gff3.gff3', 'r') as file:
	data = file.read()
gffutils.create_db(data, dbfn='chr1_db.db', force=True, from_string=True,merge_strategy="merge")

#Produces a list called "genes" which only contains the gene information section of the GFF3 file:

db = gffutils.FeatureDB('chr1_db.db', keep_order=True)
gene_features = db.features_of_type("gene")
genes = list(gene_features)

f=open("chr1_results.txt","a+")
f.write("%s\t\t\t\t%s\t%s\t\t\t%s\t%s\t%s\t\t\t%s\t\n" %("Gene_ID", "PAM_F", "gRNA_F", "GC_Per_F", "PAM_R", "gRNA_R", "GC_Per_R"))

for i in range(0,len(genes)):
	f.write("%s\n" %(genes[i].id))
#	print("%s\n" %(genes[i].id))
#	print("%s\n" %(genes[i].strand))
#If genes are located on the positive strand, the codes takes the first 300 bp of each gene:
	if genes[i].strand == "+":
		start = genes[i].start
		end = genes[i].start + 300
#Opens the chromosome1.txt file to locate each gene:
		with open("chr1.txt", "rt") as fh:
			seqs = aspairs(fh)
			for seq in seqs:
				chr_name = seq[0]
				chr_str = seq[1]
			gRNA_gene = chr_str[start:end]
#			print("+ dir: %s\n"%(gRNA_gene))
#Finds PAM sites (NGG) and gRNAs for the first 300 bp of each gene located at the positive strand:
			for i in range(23, len(gRNA_gene)):
				if gRNA_gene[i]=="G" and gRNA_gene[i-1]=="G":

					PAM_FWD = ("%s%s%s" %(gRNA_gene[i-2], gRNA_gene[i-1], gRNA_gene[i]))
					gRNA_FWD = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene[i-22],gRNA_gene[i-21],gRNA_gene[i-20],gRNA_gene[i-19],gRNA_gene[i-18],gRNA_gene[i-17],\
					gRNA_gene[i-16],gRNA_gene[i-15],gRNA_gene[i-14],gRNA_gene[i-13],gRNA_gene[i-12],gRNA_gene[i-11],gRNA_gene[i-10],gRNA_gene[i-9],gRNA_gene[i-8],gRNA_gene[i-7],\
					gRNA_gene[i-6],gRNA_gene[i-5],gRNA_gene[i-4],gRNA_gene[i-3]))
#					print("PAM_f: %s \t gRNA_f: %s" %(PAM_FWD,gRNA_FWD))
					G = 0
					C = 0
					for n in gRNA_FWD:
						if n == "G":
							G += 1
						if n == "C":
							C += 1
						else:
							continue
					GC_Per_FWD = ("%f" %(((G + C)/20)*100))
#Assessment of the uniqeness of each gRNA in the chromosome:
					if float(GC_Per_FWD) <= 50.0000 and len(re.findall(gRNA_FWD, chr_str)) == 1 and len(re.findall(r'TTT[T]+', gRNA_FWD)) < 1:
						f.write("\t\t\t\t%s\t%s\t%s\n" %(PAM_FWD, gRNA_FWD, GC_Per_FWD))

#Finds PAM and gRNA in the reverse strand:
			reversed_gRNA_gene = (''.join(reversed(gRNA_gene)))
			gRNA_gene_rev_str = []
			for i in range(0,len(reversed_gRNA_gene)):
				gRNA_gene_rev_str.append((reversed_gRNA_gene[i].translate(str.maketrans("ATCG", "TAGC"))))
			gRNA_gene_r = (''.join(gRNA_gene_rev_str))
#			print("\n- dir: %s\n" %(gRNA_gene_r))

			for j in range(23, len(gRNA_gene_r)):
				if gRNA_gene_r[j]=="G" and gRNA_gene_r[j-1]=="G":

					PAM_RVS = ("%s%s%s" %(gRNA_gene_r[j-2], gRNA_gene_r[j-1], gRNA_gene_r[j]))

					gRNA_RVS = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene_r[j-22],gRNA_gene_r[j-21],gRNA_gene_r[j-20],gRNA_gene_r[j-19],gRNA_gene_r[j-18],\
					gRNA_gene_r[j-17],gRNA_gene_r[j-16],gRNA_gene_r[j-15],gRNA_gene_r[j-14],gRNA_gene_r[j-13],gRNA_gene_r[j-12],gRNA_gene_r[j-11],gRNA_gene_r[j-10],\
					gRNA_gene_r[j-9],gRNA_gene_r[j-8],gRNA_gene_r[j-7],gRNA_gene_r[j-6],gRNA_gene_r[j-5],gRNA_gene_r[j-4],gRNA_gene_r[j-3]))
#					print("PAM_r: %s \t gRNA_r: %s" %(PAM_RVS ,gRNA_RVS))
					C = 0
					G = 0
					for n in gRNA_RVS:
						if n == "G":
							G += 1
						if n == "C":
							C += 1
						else:
							continue
					GC_Per_RVS = ("%f" %(((G + C)/20)*100))

					if float(GC_Per_RVS) <= 50.0000 and len(re.findall(gRNA_RVS, chr_str)) == 0 and len(re.findall(r'TTT[T]+', gRNA_RVS)) < 1:
						f.write("\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n" %(PAM_RVS, gRNA_RVS, GC_Per_RVS))
		continue
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if genes[i].strand == "-":
		start = genes[i].stop
		end = start - 300
		with open("chr1.txt", "rt") as fh:
			seqs = aspairs(fh)
			for seq in seqs:
				chr_name = seq[0]
				chr_str = seq[1]
			gRNA_gene = chr_str[end:start]
#			print("Forward: %s" %(gRNA_gene))
			for i in range(23, len(gRNA_gene)):
				if gRNA_gene[i] == "G" and gRNA_gene[i-1] == "G":
					PAM_FWD = ("%s%s%s" %(gRNA_gene[i-2],gRNA_gene[i-1],gRNA_gene[i]))
					gRNA_FWD = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene[i-22],gRNA_gene[i-21],gRNA_gene[i-20],gRNA_gene[i-19],gRNA_gene[i-18],\
					gRNA_gene[i-17],gRNA_gene[i-16],gRNA_gene[i-15],gRNA_gene[i-14],gRNA_gene[i-13],gRNA_gene[i-12],gRNA_gene[i-11],gRNA_gene[i-10],gRNA_gene[i-9],\
					gRNA_gene[i-8],gRNA_gene[i-7],gRNA_gene[i-6],gRNA_gene[i-5],gRNA_gene[i-4],gRNA_gene[i-3]))
#					print("-: PAM_p: %s \t gRNA_p: %s" %(PAM_RVS, gRNA_RVS))

					G = 0
					C = 0
					for n in gRNA_FWD:
						if n == "G":
							G += 1
						if n == "C":
							C += 1
						else:
							continue
						GC_Per_FWD = ("%f" %(((G + C)/20)*100))
					if float(GC_Per_FWD) <= 50.0000 and len(re.findall(gRNA_FWD, chr_str)) == 1 and len(re.findall(r'TTT[T]+', gRNA_FWD)) < 1:
						f.write("\t\t\t\t%s\t%s\t%s\n" %(PAM_FWD, gRNA_FWD, GC_Per_FWD))

#Finds the PAM and gRNAs in the reverse strand:
			reversed_gRNA_gene = (''.join(reversed(gRNA_gene)))
			gRNA_gene_rev_str = []
			for i in range(0, len(reversed_gRNA_gene)):
				gRNA_gene_rev_str.append(reversed_gRNA_gene[i].translate(str.maketrans("ATCG", "TAGC")))
			gRNA_gene_r = (''.join(gRNA_gene_rev_str))
#			print("gRNA_gene_r: %s" %(gRNA_gene_r))
			for j in range(23, len(gRNA_gene_r)):
				if gRNA_gene_r[j]=="G" and gRNA_gene_r[j-1]=="G":
					PAM_RVS = ("%s%s%s" %(gRNA_gene_r[j-2], gRNA_gene_r[j-1], gRNA_gene_r[j]))
					gRNA_RVS = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene_r[j-22],gRNA_gene_r[j-21],gRNA_gene_r[j-20],gRNA_gene_r[j-19],gRNA_gene_r[j-18],\
					gRNA_gene_r[j-17],gRNA_gene_r[j-16],gRNA_gene_r[j-15],gRNA_gene_r[j-14],gRNA_gene_r[j-13],gRNA_gene_r[j-12],gRNA_gene_r[j-11],gRNA_gene_r[j-10],\
					gRNA_gene_r[j-9],gRNA_gene_r[j-8],gRNA_gene_r[j-7],gRNA_gene_r[j-6],gRNA_gene_r[j-5],gRNA_gene_r[j-4],gRNA_gene_r[j-3]))
#					print("-: PAM_r: %s \t gRNA_r: %s" %(PAM_RVS, gRNA_RVS))

					C = 0
					G = 0
					for n in gRNA_RVS:
						if n == "G":
							G += 1
						if n =="C":
							C += 1
						else:
							continue
						GC_Per_RVS = ("%f" %(((G + C)/20)*100))
					if float(GC_Per_RVS) <= 50.0000 and len(re.findall(gRNA_RVS, chr_str)) == 0 and len(re.findall(r'TTT[T]+', gRNA_RVS)) < 1:
						f.write("\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n" %(PAM_RVS, gRNA_RVS, GC_Per_RVS))


