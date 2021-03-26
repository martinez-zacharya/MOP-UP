from __future__ import print_function, unicode_literals
import os
import pandas
import subprocess
import argparse
import multiprocessing
from BiPipelineFunctions import CutToGenome, RemoveSingletons,CodeGenomes, ExtractTitulars, ExtractFamilies, ExtractSubgroupMembers
from pyfiglet import Figlet
from time import process_time


if __name__ == "__main__":

	f = Figlet(font='speed')
	ff = Figlet(font='doh')
	print (f.renderText('MicroPipe'))


	parser = argparse.ArgumentParser()

	parser.add_argument('nameofrun',
						help='Name of the run',
						action='store')

	parser.add_argument('fasta',
						help='Enter full path to fasta file',
						action='store')

	parser.add_argument('delimiter',
						help='Delimiter between genome and gene in fasta',
						action='store')

	parser.add_argument('--miniden',
						help='Minimum percent identity to accept blast hits for building families. Default is 0.35',
						action='store',
						default= .35,
						dest='miniden')

	parser.add_argument('--minover',
						help='Minimum percent overlap to accept blast hits for building families. Default is 0.80',
						action='store',
						default= .8,
						dest='minover')

	parser.add_argument('--singleton',
						help='Add this flag to remove singletons',
						action='store_true')

	parser.add_argument('--cpu',
						help='Enter how many threads for diamond to use',
						action='store',
						default= 1,
						dest='diamondthreads')
	parser.add_argument('--iter',
						help='Enter how many infomap iterations to run',
						action='store',
						default=1000,
						dest='iters')

	parser.add_argument('--connect',
						help='Add this flag to keep only proteins that make connections in the output',
						action='store_true')

	parser.add_argument('outputfolder',
						help='Enter full path to output directory',
						action='store')

	parser.add_argument('--db',
						help = 'Type n to not use preexisting microviridae database in analysis',
						action = 'store',
						default = 'y',
						dest = 'dbcheck')


	args = parser.parse_args()

	fasta = args.fasta
	delim = args.delimiter
	minimumidentity = str(args.miniden)
	minimumoverlap = str(args.minover)
	sinless = args.singleton
	runname = args.nameofrun
	threads = str(args.diamondthreads)
	outputpath = args.outputfolder
	db = args.dbcheck
	infoiters = str(args.iters)


	if db == 'y':
		database = open('20210303.fasta', 'r')
		inputfasta = open(fasta, 'r')
		finalfasta = open('final.fasta', 'a+')
		lines1 = inputfasta.readlines()
		for line in lines1:
			finalfasta.write(line)
		lines2 = database.readlines()
		for line in lines2:
			finalfasta.write(line)

		finalfasta.close()
		database.close()
		inputfasta.close()

		final = str(os.getcwd()) + '/final.fasta'
	else:
		final = fasta

	diamondcmd = str(os.getcwd()) + '/./diamond'
	silixcmd = str(os.getcwd()) + '/silix'

	subprocess.run([diamondcmd, "makedb", "--in", final, "-d", "db"], check = True)
	subprocess.run([diamondcmd, "blastp", "-d", "db", "-q", final,"-o", "allvall.csv", "-p", 'threads'], check = True)

	outputpath = outputpath + "/output"

	with open('clusteroutput.txt', 'w') as file:
		subprocess.run(["mkdir", outputpath])
		subprocess.run([silixcmd, "-i", minimumidentity, "-r", minimumoverlap, final, "allvall.csv"], stdout = file)

	outputog = outputpath

	outputpath = outputpath + '/'

	CutToGenome('clusteroutput.txt', delim)

	if args.singleton == True:
		RemoveSingletons('CutFile.txt')
		CodeGenomes('CutFileSinless.txt')
	else:
		CodeGenomes('CutFile.txt')

	subprocess.run(['infomap', '-i', 'bipartite', '--clu', '-2', '-N', infoiters, '-s', '1', 'Coded.txt', './'], check = True)

	#Import Silix output
	clustdf = pandas.read_csv('clusteroutput.txt', delimiter='	', names=['ProteinCluster', 'Gene'])

	#Import raw InfoMap output
	df2 = pandas.read_csv('Coded.clu', delimiter = ' ', names = ["A", "B", "C"], comment = '#')

	#Import infomap input
	df1 = pandas.read_csv('Coded.txt', delimiter = ' ', names = ["A", "B"])

	# Import Silix data that has genomes to protein clusters
	if args.singleton == True:
		decodedf = pandas.read_csv('CutFileSinless.txt', delimiter = ' ', names = ["Genome","Cluster"])
	else:
		decodedf = pandas.read_csv('CutFile.txt', delimiter = ' ', names = ["Genome","Cluster"])

	decodedf = pandas.read_csv('CutFile.txt', delimiter = ' ', names = ["Genome","Cluster"])

	#Merge infomap input with genomes to protein clusters to get a table with network ID and
	#Protein clusters from infomap input, and also genomes and protein clusters from the modified
	#Silix output
	fulldf = df1.merge(decodedf, how='outer', left_index=True, right_index=True)

	#Merges the infomap output to the df with genomes, PCs and network IDs on network IDs,
	#which both dfs share. The infomap output also adds the subgroups that each genome 
	#belongs to
	df3 = fulldf.merge(df2, how='left', left_on="A", right_on="A")

	#Drops all duplicates of protein clusters, to get one genome for each subgroup,
	#which will be used for naming the subgroups
	firstuniqs = df3.drop_duplicates(subset='B_y', keep='first')
	firstuniqs = firstuniqs.drop(columns=['B_x','Cluster', 'C'])

	#This merges the new subgroup name (first genome alphabetically in the subgroup),
	#to df3, which gives each genome a subgroup name
	newdf = df3.merge(firstuniqs, how = 'left', left_on='B_y', right_on='B_y')
	df5 = newdf.drop(columns = ["A_x", "B_y", "C", "Genome_x", "B_x"])


	#for the master spreadsheet to link network id to genome to subgroup
	#This drops all duplicate entries that arise from merging
	spreadf = newdf.drop_duplicates(keep = 'first', inplace=False)
	#Drops uneccessary columns
	spreadf = spreadf.drop(columns = ["Cluster", "B_y", "B_x", "C", "A_y", "C"])
	#Drops duplicates
	spreadf = spreadf.drop_duplicates(keep='first', inplace=False)
	spreadf.to_csv(outputpath + runname + 'Master.csv', mode='w', header = ["NetworkID","Genome", "Subgroup"], index = None, sep=',')

	#For normal cytoscape visualization
	df5 = df5.drop_duplicates(keep = 'first', inplace=False)
	masterdf2 = pandas.read_csv(outputpath + runname + 'Master.csv', delimiter = ',')
	newcol2 = masterdf2['Subgroup'].value_counts().rename_axis('Unique').reset_index()
	df5 = df5.merge(newcol2, how = 'left', left_on='Genome_y', right_on='Unique')
	df5 = df5.drop(columns=['A_y', 'Unique'])
	clustdf = clustdf.drop_duplicates(subset='ProteinCluster', keep='first', inplace=False)
	df6 = df5.merge(clustdf, how='left', left_on='Cluster', right_on='ProteinCluster')
	df6 = df6.drop(columns=['ProteinCluster', 'Cluster'])
	clustdf['Annot'] = 'ProtCluster'
	clustdf = clustdf.drop(columns = ['ProteinCluster'])
	clustdf.to_csv(outputpath + runname + 'CytoscapeHelper.csv',index = None, sep=',', mode='w', header=['Protein', 'Annotation'])
	df6.to_csv(outputpath + runname + 'ForCytoscape.csv',index = None, sep=',', mode='w', header=['Subgroup', 'SubgroupCount', 'ProteinCluster'])

	subprocess.run(['rm', 'Coded.clu', 'Coded.txt', 'CutFile.txt','allvall.csv', 'db.dmnd'])
	if args.singleton == True:
		subprocess.run(['rm', 'CutFileSinless.txt'])

	path = os.path.abspath(os.getcwd())

	clustpath = str(path) + '/clusteroutput.txt'

	os.rename(clustpath, outputpath+'clusteroutput.txt')

	os.chdir(outputpath)

	subprocess.run(['mkdir', 'SubgroupMemberLists'])
	outdir = outputpath + 'SubgroupMemberLists/'
	subprocess.run(['mkdir', 'ProteinFamilies'])
	outfolder = outputpath + 'ProteinFamilies/'


	p1 = multiprocessing.Process(target=ExtractSubgroupMembers, args = (runname+'Master.csv', outdir))
	p2 = multiprocessing.Process(target=ExtractTitulars, args = (outputpath+runname+'ForCytoscape.csv', final, args.connect))
	p3 = multiprocessing.Process(target=ExtractFamilies, args = ('clusteroutput.txt', runname+'ForCytoscape.csv', final, outfolder))
		
	p1.start()
	p2.start()
	p3.start()

	p1.join()
	p2.join()
	p3.join()

	print (ff.renderText('Fin!'))
