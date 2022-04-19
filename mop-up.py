from __future__ import print_function, unicode_literals
import os
import pandas
import subprocess
import argparse
import multiprocessing
from utils import (
    CutToGenome,
    RemoveSingletons,
    CodeGenomes,
    ExtractTitulars,
    ExtractFamilies,
    ExtractSubgroupMembers,
)
from pyfiglet import Figlet
from time import process_time
from pyfaidx import Fasta


if __name__ == "__main__":

    f = Figlet(font="speed")
    ff = Figlet(font="doh")
    print(f.renderText("MOP-UP"))

    parser = argparse.ArgumentParser()

    parser.add_argument("nameofrun", help="Name of the run", action="store")

    parser.add_argument("fasta", help="Enter full path to fasta file", action="store")

    parser.add_argument(
        "delimiter", help="Delimiter between genome and gene in fasta", action="store"
    )

    parser.add_argument(
        "--miniden",
        help="Minimum percent identity to accept blast hits for building families. Default is 0.30",
        action="store",
        default=0.30,
        dest="miniden",
    )

    parser.add_argument(
        "--minover",
        help="Minimum percent overlap to accept blast hits for building families. Default is 0.80",
        action="store",
        default=0.8,
        dest="minover",
    )

    parser.add_argument(
        "--singleton", help="Add this flag to remove singletons", action="store_true"
    )

    parser.add_argument(
        "--cpu",
        help="Enter how many threads for diamond to use",
        action="store",
        default=1,
        dest="diamondthreads",
    )
    parser.add_argument(
        "--iter",
        help="Enter how many infomap iterations to run",
        action="store",
        default=1000,
        dest="iters",
    )

    parser.add_argument(
        "outputfolder", help="Enter full path to output directory", action="store"
    )

    parser.add_argument(
        "--noMicro",
        help="Add flag to not use included Microviridae database",
        action="store_true",
        dest="dbcheck",
    )
    parser.add_argument(
        "--block-size",
        help="Change Diamond Search block size to increase performance. Default is 2",
        action="store",
        default=2,
        dest="block",
    )
    parser.add_argument(
        "--sensitivity",
        help="Tune alignment sensitivity for Diamond. Choose between mid-sensitive, sensitive, more-sensitive, very-sensitive, and ultra-sensitive.",
        action="store",
        default=False,
        dest="DiamondSens",
    )

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
    blocksize = str(args.block)
    if args.DiamondSens != False:
        dmndsens = "--" + str(args.DiamondSens)
    else:
        dmndsens = args.DiamondSens

    blockarg = "-b" + blocksize

    workingDirectory = os.getcwd()

    if db == False:
        with open("final.fasta", "w") as file:
            subprocess.run(["cat", "Micro2022_02_24_input.fasta", fasta], stdout=file)
        file.close()

        final = str(os.getcwd()) + "/final.fasta"
    else:
        with open("final.fasta", "w") as finalfile:
            subprocess.run(["cat", fasta], stdout=finalfile)
        finalfile.close()
        final = str(os.getcwd()) + "/final.fasta"

    diamondcmd = str(os.getcwd()) + "/./diamond"
    silixcmd = str(os.getcwd()) + "/silix"
    genomefai = Fasta("microgenomes_20220204.fasta")

    subprocess.run([diamondcmd, "makedb", "--in", final, "-d", "db"], check=True)

    if dmndsens == False:
        subprocess.run(
            [
                diamondcmd,
                "blastp",
                "-d",
                "db",
                "-q",
                final,
                "-o",
                "allvall.csv",
                "-p",
                "threads",
                blockarg,
            ],
            check=True,
        )
    else:
        subprocess.run(
            [
                diamondcmd,
                "blastp",
                "-d",
                "db",
                "-q",
                final,
                "-o",
                "allvall.csv",
                "-p",
                "threads",
                blockarg,
                dmndsens,
            ],
            check=True,
        )

    outputpath = outputpath + "/output"

    with open("clusteroutput.txt", "w") as file:
        subprocess.run(["mkdir", outputpath])
        subprocess.run(
            [
                silixcmd,
                "-i",
                minimumidentity,
                "-r",
                minimumoverlap,
                final,
                "allvall.csv",
            ],
            stdout=file,
        )

    outputog = outputpath

    outputpath = outputpath + "/"

    CutToGenome("clusteroutput.txt", delim)

    if args.singleton == True:
        RemoveSingletons("CutFile.txt")
        CodeGenomes("CutFileSinless.txt")
    else:
        CodeGenomes("CutFile.txt")

    subprocess.run(
        [
            "infomap",
            "-i",
            "bipartite",
            "--clu",
            "-2",
            "-N",
            infoiters,
            "-s",
            "1",
            "Coded.txt",
            "./",
        ],
        check=True,
    )

    # Import Silix output
    clustdf = pandas.read_csv(
        "clusteroutput.txt", delimiter="	", names=["ProteinCluster", "Gene"]
    )

    # Import raw InfoMap output
    df2 = pandas.read_csv(
        "Coded.clu", delimiter=" ", names=["A", "B", "C"], comment="#"
    )

    # Import infomap input
    df1 = pandas.read_csv("Coded.txt", delimiter=" ", names=["A", "B"])

    # Import Silix data that has genomes to protein clusters
    if args.singleton == True:
        decodedf = pandas.read_csv(
            "CutFileSinless.txt", delimiter=" ", names=["Genome", "Cluster"]
        )
    else:
        decodedf = pandas.read_csv(
            "CutFile.txt", delimiter=" ", names=["Genome", "Cluster"]
        )

    decodedf = pandas.read_csv(
        "CutFile.txt", delimiter=" ", names=["Genome", "Cluster"]
    )

    # Merge infomap input with genomes to protein clusters to get a table with network ID and
    # protein clusters from infomap input, and also genomes and protein clusters from the modified
    # silix output
    fulldf = df1.merge(decodedf, how="outer", left_index=True, right_index=True)

    # Merges the infomap output to the df with genomes, PCs and network IDs on network IDs,
    # which both dfs share. The infomap output also adds the subgroups that each genome
    # belongs to
    df3 = fulldf.merge(df2, how="left", left_on="A", right_on="A")

    # Drops all duplicates of protein clusters, to get one genome for each subgroup,
    # which will be used for naming the subgroups
    firstuniqs = df3.drop_duplicates(subset="B_y", keep="first")
    firstuniqs = firstuniqs.drop(columns=["B_x", "Cluster", "C"])
    firstuniqs.to_csv('firstuniqs.csv')
    # This merges the new subgroup name (first genome alphabetically in the subgroup),
    # to df3, which gives each genome a subgroup name
    newdf = df3.merge(firstuniqs, how="left", left_on="B_y", right_on="B_y")
    df5 = newdf.drop(columns=["A_x", "B_y", "C", "Genome_x", "B_x"])
    df5.to_csv('df5.csv')
    # for the master spreadsheet to link network id to genome to subgroup
    # This drops all duplicate entries that arise from merging
    spreadf = newdf.drop_duplicates(keep="first", inplace=False)
    # Drops uneccessary columns
    # spreadf = spreadf.drop(columns=["B_y", "B_x", "C", "A_y", "C"])
    spreadf = spreadf.drop(columns=["Cluster", "B_y", "B_x", "C", "A_y", "C"])
    # Drops duplicates
    spreadf = spreadf.drop_duplicates(keep="first", inplace=False)
    # Imports metadata
    metadata = pandas.read_csv("Metadata.csv", comment="#", header = 1)
    metadata = metadata[["Genome designation", "Reported Source", "Source Attribution for Analyses", "GC content", "Sequence Length (nt)", "Phylum", "Class", "Order", "Family", "Genus", "Genus CRISPR prediction", "CRISPR Prediction Correct?"]]
    spreadf = spreadf.merge(metadata, how="left", left_on="Genome_x", right_on="Genome designation")
    spreadf = spreadf.drop(columns=["Genome designation"])
    spreadf.to_csv('spreadf2.csv')
    spreadf.to_csv(
        outputpath + runname + "Master.csv",
        mode="w",
        header=["NetworkID", "Genome designation", "Subgroup", "Reported Source", "Source Attribution for Analyses", "GC content", "Sequence Length (nt)", "Phylum", "Class", "Order", "Family", "Genus", "Genus CRISPR prediction", "CRISPR Prediction Correct?"],
        index=None,
        sep=",",
    )

    # For normal cytoscape visualization
    df5 = df5.drop_duplicates(keep="first", inplace=False)
    masterdf2 = pandas.read_csv(outputpath + runname + "Master.csv", delimiter=",")
    newcol2 = masterdf2["Subgroup"].value_counts().rename_axis("Unique").reset_index()
    df5 = df5.merge(newcol2, how="left", left_on="Genome_y", right_on="Unique")
    df5 = df5.drop(columns=["A_y", "Unique"])
    clustdf = clustdf.drop_duplicates(
        subset="ProteinCluster", keep="first", inplace=False
    )
    df6 = df5.merge(clustdf, how="left", left_on="Cluster", right_on="ProteinCluster")
    df6 = df6.drop(columns=["ProteinCluster", "Cluster"])
    clustdf["Annot"] = "ProtCluster"
    clustdf = clustdf.drop(columns=["ProteinCluster"])
    clustdf.to_csv(
        outputpath + runname + "CytoscapeHelper.csv",
        index=None,
        sep=",",
        mode="w",
        header=["Protein", "Annotation"],
    )
    df6.to_csv(
        outputpath + runname + "ForCytoscape.csv",
        index=None,
        sep=",",
        mode="w",
        header=["Subgroup", "SubgroupCount", "ProteinCluster"],
    )

    subprocess.run(
        ["rm", "Coded.clu", "Coded.txt", "CutFile.txt", "allvall.csv", "db.dmnd"]
    )
    if args.singleton == True:
        subprocess.run(["rm", "CutFileSinless.txt"])

    path = os.path.abspath(os.getcwd())

    clustpath = str(path) + "/clusteroutput.txt"

    os.rename(clustpath, outputpath + "clusteroutput.txt")

    os.chdir(outputpath)

    SubList = runname + "SubgroupMemberLists"
    ProtFam = runname + "ProteinFamilies"

    subprocess.run(["mkdir", SubList])
    outdir = outputpath + SubList + "/"
    subprocess.run(["mkdir", ProtFam])
    outfolder = outputpath + ProtFam + "/"
    fff = Figlet(font="poison")

    fastafai = Fasta(final)
    p1 = multiprocessing.Process(
        target=ExtractSubgroupMembers,
        args=(runname + "Master.csv", outdir, genomefai, db),
    )
    p2 = multiprocessing.Process(
        target=ExtractTitulars,
        args=(outputpath + runname + "ForCytoscape.csv", final, fastafai, runname),
    )
    p3 = multiprocessing.Process(
        target=ExtractFamilies,
        args=(
            "clusteroutput.txt",
            runname + "ForCytoscape.csv",
            final,
            outfolder,
            fastafai,
        ),
    )

    p1.start()
    p2.start()
    p3.start()

    p1.join()
    p2.join()
    p3.join()

    subprocess.run(["rm", "clusteroutput.txt"])
    if os.stat("Errors.txt").st_size == 0:
        print(ff.renderText("Done"))
    else:
        print("Pipeline Failed")

    os.chdir(workingDirectory)
    subprocess.run(
        ["rm", "final.fasta", "final.fasta.fai", "microgenomes_20220204.fasta.fai"]
    )
