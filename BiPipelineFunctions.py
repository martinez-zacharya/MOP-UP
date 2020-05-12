import os
import pandas

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
