import os
import pandas
import subprocess
import re
from Bio import SeqIO
from pyfaidx import Fasta
from tqdm import tqdm
from pyfiglet import Figlet


def CutToGenome(file, delimiter):
    f = open(file, "r")
    n = open("CutFile.txt", "a+")
    data = f.readlines()

    for line in data:
        line = line.split("	")
        line0list = line[1].split(delimiter)
        n.write(line0list[0] + " " + line[0] + "\n")

    f.close()
    n.close()


def RemoveSingletons(file):
    df = pandas.read_csv(file, delimiter=" ", names=["A", "B"])
    df = df[df.duplicated(subset=["B"], keep=False)]
    df.to_csv("CutFileSinless.txt", header=None, index=None, sep=" ", mode="a")


def CodeGenomes(file):
    n = open("Coded.txt", "a+")
    df = pandas.read_csv(file, delimiter=" ", names=["Genome", "ProteinCluster"])
    f = open(file, "r")
    maxval = df["ProteinCluster"].max()

    accdict = {}

    data = f.readlines()

    counter = maxval + 1
    for line in data:
        line = line.split(" ")
        if line[0] not in accdict:
            accdict[line[0]] = counter
            counter += 1
        n.write(str(accdict[line[0]]) + " " + line[1])
    f.close()
    n.close()


def ExtractFamilies(clustfile, cytofile, fastafile, outfolder, fastafai):
    # sequences = Fasta(fastafile)

    silix = pandas.read_csv(clustfile, delimiter="	", names=["ProteinCluster", "Gene"])

    cyto = pandas.read_csv(cytofile)

    # This is commented to keep singletons in protein family directory
    # cyto = cyto[cyto.duplicated(subset='ProteinCluster', keep=False)]

    cyto = cyto.drop(columns=["Subgroup", "SubgroupCount"])

    silixuniq = silix.drop_duplicates(
        subset="ProteinCluster", keep="first", inplace=False
    )

    cytosilix = cyto.merge(
        silixuniq, how="left", left_on="ProteinCluster", right_on="Gene"
    )

    cytosilix = cytosilix.merge(
        silix, how="left", left_on="ProteinCluster_y", right_on="ProteinCluster"
    )

    cytosilix = cytosilix.drop(columns=["Gene_x", "ProteinCluster_y"])

    cytosilix = cytosilix.drop_duplicates(keep="first", inplace=False)

    clustname = cytosilix.drop_duplicates(
        subset="ProteinCluster_x", keep="first", inplace=False
    )

    clustlist = clustname["ProteinCluster_x"].tolist()

    with (open("Errors.txt", "a+")) as error:

        try:
            for clust in tqdm(clustlist):
                newfile = open(outfolder + str(clust) + ".fasta", "a+")
                snippeddf = cytosilix[cytosilix["ProteinCluster_x"] == str(clust)]
                snippeddf = snippeddf.reset_index(drop=True)
                for i in range(len(snippeddf)):
                    newfile.write(
                        ">"
                        + str(snippeddf.loc[i, "Gene_y"])
                        + "\n"
                        + str(fastafai[str(snippeddf.loc[i, "Gene_y"])][0:])
                        + "\n"
                    )
                newfile.close()
        except Exception as e:
            error.write(str(e))
            print(e)


def ExtractTitulars(cytofile, fastafile, fastafai, nameofrun):
    df = pandas.read_csv(cytofile)

    df1 = df[df.duplicated(subset=["ProteinCluster"], keep="first")]
    proteins1 = df1["ProteinCluster"]
    proteins1 = proteins1.drop_duplicates(keep="first")
    outfasta1 = open(str(nameofrun) + "ConnectingTitularProteins.fasta", "a+")
    for protein in tqdm(proteins1):
        outfasta1.write(">" + str(protein) + "\n" + str(fastafai[protein][0:]) + "\n")
    outfasta1.close()

    proteins = df["ProteinCluster"]

    proteins = proteins.drop_duplicates(keep="first")

    outfasta = open(str(nameofrun) + "AllTitularProteins.fasta", "a+")

    for protein in tqdm(proteins):
        outfasta.write(">" + str(protein) + "\n" + str(fastafai[protein][0:]) + "\n")
    outfasta.close()


def ExtractSingletons(clustfile, fastafile):
    silix = pandas.read_csv(clustfile, delimiter="	", names=["ProteinCluster", "Gene"])

    silix = silix.drop_duplicates(subset="ProteinCluster", keep=False)

    proteins = silix["Gene"]

    outfasta = open("Singletons.fasta", "a+")

    for protein in tqdm(proteins):
        fastas = SeqIO.parse(fastafile, "fasta")
        for fasta in fastas:
            name, sequence = fasta.id, str(fasta.seq)
            if name == protein:
                SeqIO.write(fasta, outfasta, "fasta")
    outfasta.close()


def ExtractSubgroupMembers(masterfile, outfolder, genomefai, db):
    os.mkdir(outfolder + "SubgroupMemberTextFiles")
    df = pandas.read_csv(masterfile)
    grouped = df.groupby(df.Subgroup)

    subgroups = df["Subgroup"]

    subgroups = subgroups.drop_duplicates(keep="first")

    for subgroup in tqdm(subgroups):
        subgroupname = str(subgroup)
        entiresubgroup = grouped.get_group(subgroup).reset_index()
        newfile = open(
            outfolder + "SubgroupMemberTextFiles/" + subgroupname + ".txt", "a+"
        )
        for ind in entiresubgroup.index:
            newfile.write(entiresubgroup["Genome"][ind] + "\n")
        newfile.close()

    if db == False:
        os.mkdir(outfolder + "SubgroupMemberFastaFiles")
        for subgroup in tqdm(subgroups):
            subgroupname = str(subgroup)
            entiresubgroup = grouped.get_group(subgroup).reset_index()
            newfile = open(
                outfolder + "SubgroupMemberFastaFiles/" + subgroupname + ".fasta", "a+"
            )
            for ind in entiresubgroup.index:
                try:
                    newfile.write(
                        ">"
                        + str(entiresubgroup["Genome"][ind])
                        + "\n"
                        + str(genomefai[entiresubgroup["Genome"][ind]][0:])
                        + "\n"
                    )
                except KeyError:
                    pass
            newfile.close()