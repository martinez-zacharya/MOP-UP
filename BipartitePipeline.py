from __future__ import print_function, unicode_literals
import os
import pandas
import subprocess
import infomap
from BiPipelineFunctions import CutToGenome, RemoveSingletons,CodeGenomes
from pyfiglet import Figlet
from PyInquirer import prompt
from pprint import pprint


f = Figlet(font='slant')
print (f.renderText('Bipartite Pipeline'))

questions = [

	{
		'type': 'input',
		'name': 'nameofrun',
		'message': 'Please enter the name of the run'

	},
	{
		'type': 'input',
		'name': 'fasta',
		'message': 'Please enter path to file for the fasta file'

	},
	{
		'type': 'input',
		'name': 'delimiter',
		'message': 'Please enter delimiter between genome and gene in fasta file'
	},
	{
		'type': 'input',
		'name': 'miniden',
		'message': 'Please enter minimum percent identity to accept blast hits \n for building families, or press enter to use default of 35%',
		'default': '.35'
	},
	{
		'type': 'input',
		'name': 'minover',
		'message': 'Please enter minimum percent overlap to accept blast hits \n for building familes, or press enter to use default of 80%',
		'default': '.8'
	},
	{
		'type': 'confirm',
		'name': 'singleton',
		'message': 'Do you wish to remove singletons?'
	},
	{
		'type': 'input',
		'name': 'diamondthreads',
		'message': 'Please enter number of threads you want to use for Diamond'
	}
]

answers = prompt(questions)

fasta = answers['fasta']
delim = answers['delimiter']
minimumidentity = answers['miniden']
minimumoverlap = answers['minover']
sinless = answers['singleton']
runname = answers['nameofrun']
threads = answers['diamondthreads']

subprocess.run(["./diamond", "makedb", "--in", fasta, "-d", "db"])
subprocess.run(["./diamond", "blastp", "-d", "db", "-q", fasta,"-o", "allvall.csv", "-p", threads])

with open('clusteroutput.txt', 'w') as file:
	subprocess.run(["mkdir", "temp"])
	subprocess.run(["mmseqs", "createdb", fasta, "DB"])
	subprocess.run(["mmseqs", "cluster", "DB", "DB_clu", "temp"])
	subprocess.run(["mmseqs", "createtsv", "DB", "DB_clu", file])
	# subprocess.run(["silix", "-i", minimumidentity, "-r", minimumoverlap, fasta, "allvall.csv"], stdout = file)

CutToGenome('clusteroutput.txt', delim)

if answers['singleton'] == True:
	RemoveSingletons('CutFile.txt')
	CodeGenomes('CutFileSinless.txt')
else:
	CodeGenomes('CutFile.txt')

# #Change N to 100 when deployed
subprocess.run(['infomap', '-i', 'bipartite', '--clu', '-2', '-N', '100', 'Coded.txt', './'])

#Import Silix output
silixdf = pandas.read_csv('clusteroutput.txt', delimiter='	', names=['ProteinCluster', 'Gene'])

#Import raw InfoMap output
df2 = pandas.read_csv('Coded.clu', delimiter = ' ', names = ["A", "B", "C"], comment = '#')

#Import infomap input
df1 = pandas.read_csv('Coded.txt', delimiter = ' ', names = ["A", "B"])

# Import Silix data that has genomes to protein clusters
if answers['singleton'] == True:
	decodedf = pandas.read_csv('CutFileSinless.txt', delimiter = ' ', names = ["Genome","Cluster"])
else:
	decodedf = pandas.read_csv('CutFile.txt', delimiter = ' ', names = ["Genome","Cluster"])

decodedf = pandas.read_csv('CutFile.txt', delimiter = ' ', names = ["Genome","Cluster"])
runname = 'testing'

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
spreadf.to_csv(runname + 'Master.csv', mode='w', header = ["NetworkID","Genome", "Subgroup"], index = None, sep=',')


#This output is to visualize on cytoscape where subgroups only connect to PCs
#if >= 75% of the genomes in the subgroup have it
df5.to_csv("prelim50MagsHumanProTest.csv", sep=',', index=None, mode='w')
df11 = df5.drop(columns=["A_y"])
#This groups the df by genome and protein cluster and gives the counts of each PC
#in each group. Then it creates a dataframe from it
grouped = df11.groupby(['Genome_y', 'Cluster']).size()
df22 = grouped.to_frame().reset_index()
#Import Master csv
masterdf=pandas.read_csv(runname + 'Master.csv', delimiter=',')
#Gets counts for the number of genomes in each subgroup
newcol = masterdf['Subgroup'].value_counts().rename_axis('Unique').reset_index()
#Merges the counts for each subgroup with the grouped df that has the counts for
#each PC in each subgroup
df33 = df22.merge(newcol, how = 'left', left_on='Genome_y', right_on='Unique')
df33 = df33.drop(columns=['Unique'])
df33.to_csv("GroupedMagsHumanProTest.csv", sep=',', mode='w', header=["ID", "ProteinCluster", "CountofProtein", "CountOfID"], index = None)
newdf2 = pandas.read_csv('GroupedMagsHumanProTest.csv', delimiter=',')
#Takes a subset of df where subgroups only connect to PCs that 75% or more
#of their genomes have it
subsetdf = newdf2[newdf2['CountofProtein'] >= ((newdf2['CountOfID']) * 0.75)]
leftoutdf = newdf2[newdf2['CountofProtein'] < ((newdf2['CountOfID']) * 0.75)]
leftoutdf = leftoutdf.drop(columns=['CountofProtein'])
subsetdf = subsetdf.drop(columns=['CountofProtein'])
df55 = df5.drop_duplicates(keep = 'first', inplace=False)
silixdf2 = silixdf.drop_duplicates(subset='ProteinCluster', keep='first', inplace=False)
df66 = df55.merge(silixdf2, how='left', left_on='Cluster', right_on='ProteinCluster')
df66 = df66.drop(columns=['Cluster', 'A_y'])
subsetdfnew = subsetdf.merge(df66, how = 'left', left_on = 'ProteinCluster', right_on = 'ProteinCluster')
leftoutdfnew = leftoutdf.merge(df66, how='left',left_on = 'ProteinCluster', right_on = 'ProteinCluster')
subsetdfnew = subsetdfnew.drop_duplicates(subset = ['ID','ProteinCluster'], keep = 'first', inplace = False)
subsetdfnew = subsetdfnew.drop(columns=['ProteinCluster', 'Genome_y'])
leftoutdfnew = leftoutdfnew.drop_duplicates(subset = ['ID','ProteinCluster'], keep = 'first', inplace = False)
leftoutdfnew = leftoutdfnew.drop(columns=['ProteinCluster', 'Genome_y'])
subsetdfnew = subsetdfnew[['ID', 'CountOfID', 'Gene']]
leftoutdfnew = leftoutdfnew[['ID', 'CountOfID', 'Gene']]
leftoutdfnew = leftoutdfnew.drop(columns=['Gene'])
cyto75 = pandas.concat([subsetdfnew, leftoutdfnew], sort=False)
cyto75.to_csv(runname + "ForCytoscape75Percent.csv", sep=',', index=None, mode='w', header=['Subgroup', 'SubgroupCount', 'ProteinCluster'])


#For normal cytoscape visualization
df5 = df5.drop_duplicates(keep = 'first', inplace=False)
masterdf2 = pandas.read_csv(runname + 'Master.csv', delimiter = ',')
newcol2 = masterdf2['Subgroup'].value_counts().rename_axis('Unique').reset_index()
df5 = df5.merge(newcol2, how = 'left', left_on='Genome_y', right_on='Unique')
df5 = df5.drop(columns=['A_y', 'Unique'])
silixdf = silixdf.drop_duplicates(subset='ProteinCluster', keep='first', inplace=False)
df6 = df5.merge(silixdf, how='left', left_on='Cluster', right_on='ProteinCluster')
df6 = df6.drop(columns=['ProteinCluster', 'Cluster'])
silixdf['Annot'] = 'ProtCluster'
silixdf.to_csv(runname + 'CytoscapeAnnotations.csv', sep=',', index=None, mode='w')
df6.to_csv(runname + 'ForCytoscape.csv',index = None, sep=',', mode='w', header=['Subgroup', 'SubgroupCount', 'ProteinCluster'])

subprocess.run(['rm', 'prelim50MagsHumanProTest.csv', 'GroupedMagsHumanProTest.csv'])