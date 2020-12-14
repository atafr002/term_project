#!/usr/bin/env python3
from __future__ import division

import gffutils
import os, itertools, re, gzip

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

url_gff = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/746/955/GCA_001746955.1_ASM174695v1/GCA_001746955.1_ASM174695v1_genomic.gff.gz"
url_fasta = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/746/955/GCA_001746955.1_ASM174695v1/GCA_001746955.1_ASM174695v1_genomic.fna.gz"

file_gff = "GCA_001746955.1_ASM174695v1_genomic.gff.gz"
file_fasta = "GCA_001746955.1_ASM174695v1_genomic.fna.gz"

if not os.path.exists(file_gff):
	os.system("curl -O %s" %(url_gff))

if not os.path.exists(file_fasta):
	os.system("curl -O %s" %(url_fasta))

with gzip.open(file_gff, "rt") as file:
	data = file.read()

gffutils.create_db(data, dbfn='GS115_db.db', force=True, from_string=True,merge_strategy="merge")

#Defining a dictionary with chromosome IDs as dictionary keys and a list of dictionary values consisting of 
#dict.value[0] = forward chromosome strand
#dict.value[1] = reverse chromosome strand
chrs = {}
with gzip.open(file_fasta, "rt") as file:
	seqs = aspairs(file)
	for seq in seqs:
		chr_id = seq[0]
		chr_str = seq[1]
		chr_str_for = chr_str.upper()

		chr_str_r = (''.join(reversed(chr_str_for)))
		chr_str_rev = []
		for i in range(0, len(chr_str_r)):
			chr_str_rev.append(chr_str_r[i].translate(str.maketrans("ATCG", "TAGC")))
		chr_str_rev = (''.join(chr_str_rev))
		chr = []
		chr.append(chr_str_for)
		chr.append(chr_str_rev)
		chrs[str(chr_id)] = chr

#Choosing only "genes" rows from the GFF file
db = gffutils.FeatureDB('GS115_db.db', keep_order=True)
gene_features = db.features_of_type("gene")
genes = list(gene_features)

f=open("GS115_gRNA.txt","a+")
f.write("%s\t\t\t\t%s\t%s\t\t\t%s\t%s\t%s\t\t\t%s\t\n" %("Gene_ID", "PAM_F", "gRNA_F", "GC_Per_F", "PAM_R", "gRNA_R", "GC_Per_R"))

for chr_id in chrs.keys():
	chr_str_for = chrs[str(chr_id)][0]
	chr_str_rev = chrs[str(chr_id)][1]
	print("%s\n" %(chr_id))

	for i in range(0,len(genes)):
		if genes[i].seqid == chr_id:
			f.write("%s\n%s\n" %(genes[i].id, chr_id))


#If genes are located on the positive strand, the codes takes the first 300 bp of each gene:
			if genes[i].strand == "+":
				start = genes[i].start
				end = genes[i].start + 300

				gRNA_gene = chr_str_for[start:end]

#Finds PAM sites (NGG) and gRNAs for the first 300 bp of each gene located at the positive strand:
				for i in range(23, len(gRNA_gene)):
					if gRNA_gene[i]=="G" and gRNA_gene[i-1]=="G":

						PAM_FWD = ("%s%s%s" %(gRNA_gene[i-2], gRNA_gene[i-1], gRNA_gene[i]))
						gRNA_FWD = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene[i-22],gRNA_gene[i-21],gRNA_gene[i-20],gRNA_gene[i-19],gRNA_gene[i-18],gRNA_gene[i-17],\
						gRNA_gene[i-16],gRNA_gene[i-15],gRNA_gene[i-14],gRNA_gene[i-13],gRNA_gene[i-12],gRNA_gene[i-11],gRNA_gene[i-10],gRNA_gene[i-9],gRNA_gene[i-8],gRNA_gene[i-7],\
						gRNA_gene[i-6],gRNA_gene[i-5],gRNA_gene[i-4],gRNA_gene[i-3]))


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
						if float(GC_Per_FWD) <= 50.0000 and len(re.findall(r'TTT[T]+', gRNA_FWD)) < 1 and len(re.findall(gRNA_FWD, chr_str_for)) == 1 and\
						len(re.findall(gRNA_FWD, chr_str_rev)) == 0:
							f.write("\t\t\t\t%s\t%s\t%s\n" %(PAM_FWD, gRNA_FWD, GC_Per_FWD))

#Finds PAM and gRNA in the reverse strand:
				reversed_gRNA_gene = (''.join(reversed(gRNA_gene)))
				gRNA_gene_rev = []
				for i in range(0,len(reversed_gRNA_gene)):
					gRNA_gene_rev.append((reversed_gRNA_gene[i].translate(str.maketrans("ATCG", "TAGC"))))
				gRNA_gene_r = (''.join(gRNA_gene_rev))



				for j in range(23, len(gRNA_gene_r)):
					if gRNA_gene_r[j]=="G" and gRNA_gene_r[j-1]=="G":

						PAM_RVS = ("%s%s%s" %(gRNA_gene_r[j-2], gRNA_gene_r[j-1], gRNA_gene_r[j]))

						gRNA_RVS = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene_r[j-22],gRNA_gene_r[j-21],gRNA_gene_r[j-20],gRNA_gene_r[j-19],gRNA_gene_r[j-18],\
						gRNA_gene_r[j-17],gRNA_gene_r[j-16],gRNA_gene_r[j-15],gRNA_gene_r[j-14],gRNA_gene_r[j-13],gRNA_gene_r[j-12],gRNA_gene_r[j-11],gRNA_gene_r[j-10],\
						gRNA_gene_r[j-9],gRNA_gene_r[j-8],gRNA_gene_r[j-7],gRNA_gene_r[j-6],gRNA_gene_r[j-5],gRNA_gene_r[j-4],gRNA_gene_r[j-3]))

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

						if float(GC_Per_RVS) <= 50.0000 and len(re.findall(r'TTT[T]+', gRNA_RVS)) < 1 and len(re.findall(gRNA_RVS, chr_str_rev)) == 1 and\
						len(re.findall(gRNA_RVS, chr_str_for)) == 0:
							f.write("\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n" %(PAM_RVS, gRNA_RVS, GC_Per_RVS))

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			if genes[i].strand == "-":
				start = genes[i].stop
				end = start - 300

				gRNA_gene = chr_str_for[end:start]

				for i in range(23, len(gRNA_gene)):
					if gRNA_gene[i] == "G" and gRNA_gene[i-1] == "G":
						PAM_FWD = ("%s%s%s" %(gRNA_gene[i-2],gRNA_gene[i-1],gRNA_gene[i]))
						gRNA_FWD = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene[i-22],gRNA_gene[i-21],gRNA_gene[i-20],gRNA_gene[i-19],gRNA_gene[i-18],\
						gRNA_gene[i-17],gRNA_gene[i-16],gRNA_gene[i-15],gRNA_gene[i-14],gRNA_gene[i-13],gRNA_gene[i-12],gRNA_gene[i-11],gRNA_gene[i-10],gRNA_gene[i-9],\
						gRNA_gene[i-8],gRNA_gene[i-7],gRNA_gene[i-6],gRNA_gene[i-5],gRNA_gene[i-4],gRNA_gene[i-3]))

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

						if float(GC_Per_FWD) <= 50.0000 and len(re.findall(r'TTT[T]+', gRNA_FWD)) < 1 and len(re.findall(gRNA_FWD, chr_str_for)) == 1 and\
						len(re.findall(gRNA_FWD, chr_str_rev)) == 0:
							f.write("\t\t\t\t%s\t%s\t%s\n" %(PAM_FWD, gRNA_FWD, GC_Per_FWD))

#Finds the PAM and gRNAs in the reverse strand:
				reversed_gRNA_gene = (''.join(reversed(gRNA_gene)))
				gRNA_gene_rev_str = []
				for i in range(0, len(reversed_gRNA_gene)):
					gRNA_gene_rev_str.append(reversed_gRNA_gene[i].translate(str.maketrans("ATCG", "TAGC")))
				gRNA_gene_r = (''.join(gRNA_gene_rev_str))

				for j in range(23, len(gRNA_gene_r)):
					if gRNA_gene_r[j]=="G" and gRNA_gene_r[j-1]=="G":

						PAM_RVS = ("%s%s%s" %(gRNA_gene_r[j-2], gRNA_gene_r[j-1], gRNA_gene_r[j]))

						gRNA_RVS = ("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s" %(gRNA_gene_r[j-22],gRNA_gene_r[j-21],gRNA_gene_r[j-20],gRNA_gene_r[j-19],gRNA_gene_r[j-18],\
						gRNA_gene_r[j-17],gRNA_gene_r[j-16],gRNA_gene_r[j-15],gRNA_gene_r[j-14],gRNA_gene_r[j-13],gRNA_gene_r[j-12],gRNA_gene_r[j-11],gRNA_gene_r[j-10],\
						gRNA_gene_r[j-9],gRNA_gene_r[j-8],gRNA_gene_r[j-7],gRNA_gene_r[j-6],gRNA_gene_r[j-5],gRNA_gene_r[j-4],gRNA_gene_r[j-3]))

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

						if float(GC_Per_RVS) <= 50.0000 and len(re.findall(r'TTT[T]+', gRNA_RVS)) < 1 and len(re.findall(gRNA_RVS, chr_str_rev)) == 1  and\
						len(re.findall(gRNA_RVS, chr_str_for)) == 0:
							f.write("\t\t\t\t\t\t\t\t\t\t%s\t%s\t%s\n" %(PAM_RVS, gRNA_RVS, GC_Per_RVS))


