Version : 0.9.5


**SETUP**

1. This has only been tested on Ubuntu, but it might work on other Linux distributions.

2. Input ```$ git clone https://github.com/martinez-zacharya/MOP-UP``` into the command line

3. You need to install various things including

  	1. If using Ubuntu, build-essential ```$ sudo apt-get install build-essential```. You might already have this installed.
  	2. SiLiX (https://lbbe-web.univ-lyon1.fr/fr/SiLix)
		1. Make sure to place the silix file into the same directory as mop-up.py and make it exectuable (```$ chmod +x silix```)
	3. Miniconda (https://docs.conda.io/en/latest/miniconda.html)
4. Enter the MOP-UP directory, enter ```$conda activate```
5. Then you can enter ```$ conda env create -f environment.yml``` to create the necessary conda environment


6. You need a protein fasta file with 
	1. No spaces in the headers
	2. A delimiter that serves as a cutoff between genome name and gene name, where to the left of the delimiter is the genome and to the right is the gene. In the example below, the delimiter is "@". Also, make sure the delimiter is a valid character that can be used in a filename, unlike "&". We recommend @.

	> \>ClusterI_0_SCPE01000001_Alistipes_sp@gene_000001
	> MEKLKALLTSKKFWTLVAAIVAALTAFFTTSCTGYLKFRREGVHHDTVRYEQVIKHKNYSAWLSNQIDRSSWRRPMMLSVCSSGTAFLSRSSGVSISSYLISPPLIHSVLPSILSSSINSVSFITLNRNCLEQITPLFTSNLPKRKSLGVCLALMMIGVPLSRVSPSSLDERRRGGRKGRPRPKIRNIPIGGSHL
	
	3. No duplicated gene names



**RUNNING**

To run, first make sure that diamond, silix and the python scripts that come with the package are in your curret working directory. The required parameters for the pipeline are

	1. Name of the run
	2. The full path to your protein fasta file
	3. The delimiter that separates the cutoff between Genome and Gene, such as '@' (@ is the default and we recommend using that)
	4. The full path to your desired output directory
	
There are several optional arguments you can add, see below. The most important is --miniden, which at 30 is used to find clusters of microviral families, and at 50 is used to define genera. You should generally run MOP-UP for both. 

	1. --miniden Minimum percent identity to accept blast hits for building families. Default is 35% # PAUL NOTE MAKE THIS 30%!
	2. --minover Minimum percent overlap to accept blast hits for building families. Default is 80%
	3. --cpu Number of threads you want Diamond to use. Default is all available threads
	4. --noMicro By adding this argument, you can elect to leave out the Microviridae database from the analysis. MOP-UP will then only work with the user-provided files
	5. --singleton Add this flag to remove singletons
	6. --iter How many infomap iterations to run. Default number is 1000. Fewer runs result in more 
	7. --block-size Change block size for Diamond to increase performance depending on available RAM
	8. --sensitivity Tune alignment sensitivity for Diamond. Options include mid-sensitive, sensitive, more-sensitive, very-sensitive, and ultra-sensitive

An example command

```$ python3 mop-up.py TestRun /path/to/Bipartite_pipeline.py/Input.fasta @ /path/to/TestRunOutput --cpu 12```

**OUTPUT**

Here is a brief description of the MOP-UP output. nameofrunForCytoscape.csv and nameofrunMaster.csv are the most important.

	1. nameofrunAllTitularProteins.fasta - One 'representative' protein from each protein cluster in a single fasta file. Used to reduce redundancy in  
					      identification via BLAST and such
	2. nameofrunConnectingTitularProteins.fasta - One 'representative' protein from each protein cluster that connects to another cluster through a subgroup in a 
	  				      single fasta file. Used to quickly identify proteins shared betweeen more distantly related phages.
	3. nameofrunForCytoscape.csv - A network file for viewing with a program such as cytoscape
	4. nameofrunCytoscapeHelper.csv - When added to Cytoscape with the network file, allows you to differentiate between the protein cluster and subgroup nodes
	5. nameofrunMaster.csv - The master file that lists all of the genomes in the Micropipe run and their respective network ID and subgroup (i.e. group of   
	                             closely related genomes)
	6. nameofrunProteinFamilies - A directory filled with fasta files for each protein cluster. Useful for alignments.
	7. nameofrunSubgroupMembers - A directory that contains text files for each subgroup with the subgroup members listed. Can be used to extract groups of closely related phages from 

**HOW TO VIEW OUTPUT IN CYTOSCAPE**

	1. Download Cytoscape from https://cytoscape.org/ and install it
	2. Open Cytoscape, File > Import Network From File, open nameofrunForCytoscape
	3. Click on the subgroup column and change it to Source Node (green dot), Protein Cluster column to Target Node (red target), click ok and wait for the 
	   network to load
	4. Once it is loaded, you want to make this more viewable by removing all proteins that are just found exclusively in a single phage group - they clutter up 
	   the graph. Click the Filter tab on the left side of the screen, add new condition, degree filter, and set it from 0 to 1 inclusive. Then click Edit at the            top bar and select "Delete Selected Nodes and Edges". Rebuild the network layout by going to the Layout tab at the top, clicking Prefuse Force Directed CL 
	   Layout > none (or whatever layout you prefer.
	5. Now a network should be created and you can play around with it. It is advisable to go to Filter > Column Filter > name and select all the proteins (they  
	   should have an @ in their names by default). With these selected go to the Style Tab on the left and change their colour/shape. You can also do the same 	       with VP1 proteins specifically, which tend to be at the center of large clusters. 

**HOW TO INTERPRET OUTPUT**

Check the nameofrunMaster.csv file to see what "subgroup" your genomes are falling into. If you ran MOP-UP at 50% miniden, they might be in a subgroup of 	        closely related genomes of syntenic gene order (operationally we call those genera, if you ran it at 30% miniden this group will include more distantly 	   related genomes too). Either way your genomes might be grouped with taxonomically defined phages and you can determine what yours are through the names of their close relatives. This could be something like FamilyX_subgroup_name-of-genome, for example Family3_Shukshmavirinae_MG945270 for a member of the Shukshmavirinae subfamily in Family 3. If there are no taxonomically assigned phages in the subgroup, search for your genomes in the Cytoscape network. Most likely, they are connected to a family-defining VP1, if not at the minid =50 network then at the minid=30 network. You can figure that out visually, or by right clicking on your subgroup of choice and clicking Select > Select first neighbours. It will highligh the proteins it is connected to, leading you to the VP1 of interest. If your genome is a singleton unconnected to any VP1s, you will have to resort to phylogenetics again. The file nameofrunConnectingTitularProteins.fasta might be useful since it will contain all required family-defining and singleton VP1s for tree building (among others). Note that many singleton genomes still clearly fall into families and are not automatically new families or whatnot. Instances of that are indicated as FamilyX_associated in the MOP-UP database. 
	   
