#!/usr/bin/env python

def MC1R_genes(position):
	"""	MC1R gene
	melanocortin 1 receptor
	https://ghr.nlm.nih.gov/gene/MC1R#location
	Molecular Location: base pairs 89,917,879 to 89,920,977 on chromosome 16
	(Homo sapiens Annotation Release 108, GRCh38.p7)"""

	position = int(position)
	# if position >= 89917879 and position <= 89920977:
	if position >= 89900000 and position <= 89930000:
		return True


dna_path = "/Users/rgroussman/Desktop/ryan/aphid/genome/chr_16.txt"
dna = open(dna_path, 'r')

"""rsid    chromosome      position        allele1 allele2
rs190214723     1       693625  T       T
rs3131972       1       752721  A       G
rs12562034      1       768448  G       G
"""

# grep -E "\t16\t" AncestryDNA.txt > chr_16.txt

for line in dna:
	if line.startswith('#') == True:
		pass
	else:
		line_elts = line.split("\t")
		rsid = line_elts[0]
		chromosome = line_elts[1]
		position = line_elts[2]
		allele1 = line_elts[3]
		allele2 = line_elts[4].strip()
		if MC1R_genes(position) == True:
			print line




# get redhair genes:
