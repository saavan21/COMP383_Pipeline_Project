# COMP383_Pipeline_Project

**Description/Task:**
We are interested in looking at the Human herpesvirus 5 transcriptomes and build a pipeline to analyze the HCMV. 

**Steps Taken (track 2):**

1. Retrieve the following transcriptomes from two patient donors from SRA and convert to paired-end fastq files. Donor 1 (2dpi), Donor 1 (6dpi), Donor 3 (2dpi), Donor 3 (6dpi). Use wget to retrieve all files 
2. To answer the question: Which strains are most similar to these patient samples, we will need to assemble transcriptome reads. Using Bowtie 2, we will create a index for HCMV (NCBI accession NC_006273.2). Only save the reads that map to the HCMV index for use in the assembly. The PipelineProject.log file will have the number of reads in each transcriptome before and after the Bowtie2 mapping. 
3. The Bowtie2 output will be used to assemble all four transcriptomes together to produce one assembly via SPAdes. The PipelineProject.log file will have the SPAdes command that was used. 
4. Create Python code to calculate the number of contigs that are > 1000 in length and store the output in the PipelineProject.log file. Additionally, create python code that will calculate the length of the assembly and write the number of bp in the assembly in the PipelineProject.log.
5. Create python code that will retrieve the longest contain from the SPAdes assembly and use blast+ input to query the nr nucleotide database. The blast+ should only keep the best alignment for the sequences. The results should be in a table format in the PipelineProject.log file. 

Requirements to run the Pipeline: 

Linux/Unix, Biopython, Python3, SRA Kit (perform fastq-dump), SPAdes, Bowtie2, Blast+

Run the Pipeline: 
To run the script, first clone the project from GitHub (https://github.com/saavan21/COMP383_Pipeline_Project.git) and type the command to run the PipelineProject_main.py which will run through the full pipeline and store data in the PipelineProject.log file. Make sure that the following dependencies are installed on the terminal or server. There are test data that can be used to run the pipeline faster. The sequence.fasta file will be needed to run the blast function. 

- Import os
- Import logging 
- from Bio import SeqIO
- from Bio import Entrez
- from Bio import SearchIO
- from Bio.Seq import Seq 
