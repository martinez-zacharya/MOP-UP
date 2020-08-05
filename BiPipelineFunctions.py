import os
import pandas
from Bio import SeqIO

def CutToGenome(file, delimiter):
	f = open(file, "r")
	n = open('CutFile.txt', 'a+')
	data = f.readlines()

	for line in data:
		line = line.split('	')
		line0list = line[1].split(delimiter)
		n.write(line0list[0] + ' ' + line[0] + '\n')

	f.close()
	n.close()


def RemoveSingletons(file):
	df = pandas.read_csv(file, delimiter = ' ', names = ["A", "B"])
	df = df[df.duplicated(subset=["B"], keep = False)]
	df.to_csv('CutFileSinless.txt', header = None, index = None, sep=' ', mode='a')


def CodeGenomes(file):
	n = open("Coded.txt", "a+")
	df = pandas.read_csv(file, delimiter = ' ', names = ['Genome', 'ProteinCluster'])
	f = open(file, "r")
	maxval = df['ProteinCluster'].max()

	accdict = {}

	data = f.readlines()

	counter = maxval + 1
	for line in data:
		line = line.split(' ')
		if line[0] not in accdict:
			accdict[line[0]] = counter
			counter += 1
		n.write(str(accdict[line[0]]) + ' ' + line[1])
	f.close()
	n.close()

def ExtractFamilies(clustfile, cytofile, fastafile, outfolder):
	silix = pandas.read_csv(clustfile, delimiter = '	', names = ['ProteinCluster', 'Gene'])

	cyto = pandas.read_csv(cytofile)

	cyto = cyto[cyto.duplicated(subset='ProteinCluster', keep=False)]

	cyto = cyto.drop(columns=['Subgroup', 'SubgroupCount'])

	silixuniq = silix.drop_duplicates(subset='ProteinCluster', keep='first', inplace=False)

	cytosilix = cyto.merge(silixuniq, how = 'left', left_on = 'ProteinCluster', right_on = 'Gene')

	cytosilix = cytosilix.merge(silix, how = 'left', left_on = 'ProteinCluster_y', right_on = 'ProteinCluster')

	cytosilix = cytosilix.drop(columns = ['Gene_x', 'ProteinCluster_y'])

	cytosilix = cytosilix.drop_duplicates(keep = 'first', inplace=False)

	clustname = cytosilix.drop_duplicates(subset = 'ProteinCluster_x', keep='first', inplace=False)

	clustlist = clustname['ProteinCluster_x'].tolist()


	for clust in clustlist:
		newfile = open(outfolder + str(clust)+'.fasta', 'a+')
		snippeddf = cytosilix[cytosilix['ProteinCluster_x'] == str(clust)]
		snippeddf = snippeddf.reset_index(drop=True)
		for i in range(len(snippeddf)):
			fasta_sequences = SeqIO.parse(fastafile, "fasta")
			for fasta in fasta_sequences:
				name, sequence = fasta.id, str(fasta.seq)
				if str(clust) == snippeddf.loc[i, "ProteinCluster_x"] and name == snippeddf.loc[i, "Gene_y"]:
					SeqIO.write(fasta, newfile, "fasta")
		newfile.close()

def ExtractTitulars(cytofile, fastafile, yn):
	df = pandas.read_csv(cytofile)

	#y = yes, only keep proteins that make connections
	#n = no, keep all proteins
	if yn == True:
		df1 = df[df.duplicated(subset=['ProteinCluster'], keep='first')]
		proteins1 = df1['ProteinCluster']
		proteins1 = proteins1.drop_duplicates(keep='first')
		outfasta1 = open('ConnectingTitularProteins.fasta', 'a+')
		for protein in proteins1:
			fastas = SeqIO.parse(fastafile, "fasta")
			for fasta in fastas:
				name, sequence = fasta.id, str(fasta.seq)
				if name == protein:
					SeqIO.write(fasta, outfasta, "fasta")
		outfasta1.close()
	

	proteins = df['ProteinCluster']

	proteins = proteins.drop_duplicates(keep='first')

	outfasta = open('AllTitularProteins.fasta', 'a+')

	for protein in proteins:
		fastas = SeqIO.parse(fastafile, "fasta")
		for fasta in fastas:
			name, sequence = fasta.id, str(fasta.seq)
			if name == protein:
				SeqIO.write(fasta, outfasta, "fasta")
	outfasta.close()

def ExtractSingletons(clustfile, fastafile):
	silix = pandas.read_csv(clustfile, delimiter = '	', names = ['ProteinCluster', 'Gene'])

	silix = silix.drop_duplicates(subset='ProteinCluster', keep = False)

	proteins = silix['Gene']

	outfasta = open('Singletons.fasta', 'a+')

	for protein in proteins:
		fastas = SeqIO.parse(fastafile, "fasta")
		for fasta in fastas:
			name, sequence = fasta.id, str(fasta.seq)
			if name == protein:
				SeqIO.write(fasta, outfasta, "fasta")
	outfasta.close()

