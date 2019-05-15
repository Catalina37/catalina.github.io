# catalina.github.io
 ---------------------------------------------------------------------------------------------------------------------------
#Barcount: written by Stephan Kamrad (stephan.kamrad.15@ucl.ac.uk).
Python package for finding barcodes within sequencing reads and matching to reference database.
barcount_test.py is a python package containing in-built files to test barcount. 

#BarSeqTools: written by StJohn Townsend (stjohn.townsend.11@ucl.ac.uk).
R package containing utility functions for barcode matching and visualisation.

#Bioneer V5.0 characterisation: folder where all characterisation related files are found. 
 #Browser: folder containing all the required files to visualise the genes in the Browser.
 
 ---------------------------------------------------------------------------------------------------------------------------

# Instructions for installing & running BarSeqTools & Barcount
BarSeqTools is a custom-built R package containing utility function to characterise the 
Bioneer version 5.0 of the fission yeast deletion collection library by barcode matching 
and visualisation.

### Install package:BarSeqTools/Barcount
------------------------------------------------------------------------------------------
1. Open R studio
2. Install devtools package
   -install.packages("devtools")
3. Load devtools package
   -library(devtools)
4. Install github package
   -install_github("Catalina37/catalina.github.io", subdir = "package name ")
------------------------------------------------------------------------------------------
The following packages will need to be available. You can install them as follows:

------------------------------------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")  <br/>
BiocInstaller::biocLite("BiocGenerics") <br/>
BiocInstaller::biocLite("Biostrings") <br/>
BiocInstaller::biocLite("DECIPHER") <br/>
BiocInstaller::biocLite("GenomicFeatures") <br/>
BiocInstaller::biocLite("GenomicRanges") <br/>
BiocInstaller::biocLite("Gviz") <br/>
BiocInstaller::biocLite("IRanges") <br/>
BiocInstaller::biocLite("S4VectorsGen") <br/>

conda install -c bioconda bedops <br/>
conda install -c conda-forge sed  <br/>

Please, note for sed usage on Mac you need to change the '-i' flag. 
This command is run in MAGIC but can be edited if required using fix(MAGIC).
This opens a window that allowing you to add "" following the '-I' flag, like so:
sed -i "" 's/\(@.*2:N:0:\).*/\1/' FileName


### Check the following R package dependencies are available:
------------------------------------------------------------------------------------------
-shiny <br/>
-Gviz <br/>
-S4Vectors <br/>
-tats4 <br/>
-BiocGenerics <br/>
-parallel <br/>
-IRanges <br/>
-GenomicRanges <br/>
-GenomeInfoDb <br/>
-grid <br/>
-GenomicFeatures <br/>
-AnnotationDbi <br/>
-Biobase <br/>

------------------------------------------------------------------------------------------
BarSeqTools is ready to use!

------------------------------------------------------------------------------------------
### To perform the deletion collection characterisation, check the following dependencies 
### are available & installed on your machine:
------------------------------------------------------------------------------------------
1. Anaconda <br/>
For details on how to install it, check out this link: https://docs.anaconda.com/anaconda
/install/

2. Barcount <br/>
A custom-built Python script used to search for the barcodes within the read using flanks 
on either side of the barcode. Barcount was developed in python2.

-----------------
Follow the instructions below to install Barcount:
- download/clone package from Github <br/>
- ensure python 2 is set as your default, otherwise in your terminal type: <br/>
  alias python=python2 <br/>
  python -V <br/>
- install pandas dependency <br/>
  python -m pip install pandas <br/>

- install biopython dependency  <br/>
  python -m pip install biopython <br/>

- install python-Levenshtein package dependency <br/>
  python -m pip install python-Levenshtein <br/>

- add Barcount to your path <br/>
-----------------

3. Bedtools <br/>
Install with conda: <br/>
install -c bioconda bedtools

4. Bowtie2 <br/>
Install with conda: <br/>
conda install -c bioconda bowtie2

5. fastqCombinedPairedEnd.py <br/>
Python script used to filter read 2 based on the header information of read 1. 
Read 1 header contains the barcode, thus matching read 1 barcode containing reads to 
their respective read 2 reads. The script can be downloaded from Github code page.

6. Fastq-grep <br/>
The package is part of fastx-toolkit. 
Install with conda:  <br/>
conda install -c bioconda fastq-tools <br/>

7. Fastx_trimmer  <br/>
The package is part of fastq-tools.   <br/>
Install with conda:  <br/>
conda install -c biobuilds fastx-toolkit <br/>

8. Samtools <br/>
Install with conda:  <br/>
conda install -c bioconda samtools  <br/>

9. PEAR  <br/>
For details on how to install it, check out this link: 
https://cme.h-its.org/exelixis/web/software/pear/doc.html)

10. Python <br/>
For details on how to install it, check out this link: https://www.python.org/downloads/
------------------------------------------------------------------------------------------


### In addition to the above, the following file dependencies are also required:
------------------------------------------------------------------------------------------

1. All protein coding genes in gtf file format.

2. A bedfile format of the Bioneer library containing all the genomic coordinates and 
strands for all of the protein coding genes.

3. An empty csv file required to run barcount since barcount was made to run with the 
barcode_table argument. 
Note, the file is not being used thus, the overall result is not affected by this.

4. The indexed S.pome genome folder.
------------------------------------------------------------------------------------------
