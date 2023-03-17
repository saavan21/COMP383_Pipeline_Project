#import packages needed for the code
import os
import logging
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW

#set varibale directory to current directory
directory = os.getcwd()
#changes current directory to sepcified directory
os.chdir(directory)

#creates and opnes a log file and will add all outputs as code runs
log_file = open('PipelineProject.log', 'a')

#create a list 'sequences' that will store all accession numbers for the files to be retrieved
sequences = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

#define function 'getSeq' to grab the sequences needed
def getSeq(seq_url):
    #wget command to download the sequence given a URL
    grabSeq = 'wget ' + seq_url
    #fast-dump command to uncompress and split the sequence files
    uncomp = 'fastq-dump -I --split-files ' + seq_url.split('/')[-1]
    #os command to download the files
    os.system(grabSeq)
    #os command to uncompress and split the sequnce files using fastq-dump
    os.system(uncomp)

#list of URL's to download the sequences
seq_urls = [
    'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030',
    'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033',
    'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044',
    'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045'
]
#for loop to go through sequence URl's and download the sequence data for each url
for seq_url in seq_urls:
    getSeq(seq_url)

#define the function testData to create the test data
def testData(file):
    #Get the current working directory
    current = os.getcwd()

    #for loop to go through each filename in the input file list
    for filename in file:
        #open the input file for reading 
        fileInput = open(current+'/'+filename,'r')
        #open the output file that will write output and add testdata_ at the front
        fileOutput = open(current+'/testdata_' + filename, 'w')

        #Read the input file and plit into a list of lines
        fileFormat = fileInput.read().strip().split('\n')

        #loop through the first 40,0000 ines of the input data 
        for i in range(40000):
            #output each line to the output file
            fileOutput.write(fileFormat[i]+'\n')

#List of file names that testData will take and process
newTestData = ['SRR5660030_1.fastq','SRR5660030_2.fastq', 'SRR5660033_1.fastq', 'SRR5660033_2.fastq', 'SRR5660044_1.fastq', 'SRR5660044_2.fastq','SRR5660045_1.fastq','SRR5660045_2.fastq']

testData(newTestData)

#Import the Entrez module from Biopython
def transcriptome():
    #Set email for Entrez
    Entrez.email = 'spatel92@luc.edu'

    #Open files for writing the output
    Fastaout = open("NC_006273.fasta",'w')
    CDSout = open("NC_006273_CDS.fasta",'w')

    #this will fetch the FASTA records for NC_006273
    title = Entrez.efetch(db='nucleotide',id='NC_006273', rettype='fasta')
    fasta_record = list(SeqIO.parse(title,'fasta'))

    #this will write the fasta record to a file
    Fastaout.write('>' + str(fasta_record[0].description) + '\n' +str(fasta_record[0].seq))
    Fastaout.close()

    #this will fetch all the genbank records for NC_006273
    GenBank_record = Entrez.efetch(db='nucleotide', id='NC_006273', rettype='gb', retmode='text')

    #count the number of CDS in the genbank record put into a file
    count = 0
    #for loop to go over reacrods in the genbank file
    for rec in SeqIO.parse(GenBank_record, 'genbank'):
        #go over each feature listed in the record 
        for i in rec.features:
            #check if the feature is a CDS
            if i.type == 'CDS':
                #if it is then count the CDS found
                count = count + 1
                #this will write the protein ID and the CDS sequence to the CDS output file 
                CDSout.write('>' + str(i.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'", "") + '\n' + str(i.location.extract(rec).seq) + '\n')
    #close the output file
    CDSout.close()
    #print the count 
    print(count)

#call function to execute the code
result = transcriptome()

#log_file.write('The HCMV genome NC_006273) has ' + str(result)+ 'CDS ' + '\n')
#print('THE HCMV genome (NC_006273) has ' + str(result)+ ' CDS. ')

#define a function that runs the bowtie2
def bowtie2(sequence):
    #build the bowtie2 index using the NC_006273 reference sequence 
    bowtie2_input = 'bowtie2-build ./NC_006273.fasta NC_006273_index'
    #call os system for file above
    os.system(bowtie2_input)

    #this will be the input sequence using the NC_006273 index
    bowtie2_index = 'bowtie2 -x NC_006273_index -1 '+'testdata_' +sequence+ '_1.fastq -2 '+ 'testdata_'+ sequence +'_2.fastq -S tmp.sam --al-conc'+ ' testdata_'+sequence+ '_bowtie2_%.fastq'
    os.system(bowtie2_index)

#for loop to go over each sequence in the sequence list and run Bowtie2 on each one
for i in sequences:
    bowtie2(i)

#define mapReads function 
def mapReads(sequence):
    #Define empty string 
    title = ''
    #if/elif statements to assign a string value to the title based on the value of the sequence
    if sequence == 'SRR5660030':
        title = 'Donor 1 (2dpi)'
    elif sequence == 'SRR5660033':
        title = 'Donor 1 (6dpi)'
    elif sequence == 'SRR5660044':
        title = 'Donor 3 (2dpi)'
    else: 
        title = 'Donor 3 (6dpi)'

    #open both files and assign to a variable 
    original_File1 = open('testdata_' + sequence + '_1.fastq')
    original_File2 = open('testdata_' + sequence + '_2.fastq')

    #define two variables with inital count at 0
    intial= 0 
    count = 0

    #Use a for loop to go over the lines in the originalfile1 and increment by 1 for each line
    for l in original_File1:
        intial = intial + 1
    #Use a for loop to go over the lines in the originalfile2 and increment by 1 for each line 
    for l in original_File2:
        count = count + 1

    #calculate the total number of read pairs before Bowtie2 filtering and assign it to countBefore
    countBefore = (intial+count)/8
    resultFile1 = open('testdata_' + sequence + '_bowtie2_1.fastq')
    resultFile2 = open('testdata_' + sequence + '_bowtie2_2.fastq')

    #Define two variables called with an initial value of 0
    count_first = 0
    count_second = 0

    #for loop to go over lines in the resultFile1 and add by 1
    for l in resultFile1:
        count_first = count_first + 1

    #for loop to go over lines in the resultFile2 and add by 1
    for l in resultFile2:
        count_second = count_second + 1

    #calculate the total number of read pairs after Bowtie2 filtering and assign it to count_final
    count_final = (count_first + count_second)/8

    #print output in the log file formated
    log_file.write(str(title)+' had ' +str(countBefore) + ' read pairs before Bowtie2 filtering and '+ str(count_final)+ ' pairs after.'+'\n')

for i in sequences:
    mapReads(i)

#define a function to run Spades
def SPAdes(files):
    #assign the four outputs file names awith variables 
    file1 = files[0]
    file2 = files[1]
    file3 = files[2]
    file4 = files[3]

    #command that will run spades 
    runSPAdes = 'spades.py -k 55,77,99,127 -t 2 --only-assembler --pe1-1 '+'testdata_'+file1+'_bowtie2_1.fastq --pe1-2 '+'testdata_'+file1+'_bowtie2_2.fastq --pe2-1 '+'testdata_'+file2+'_bowtie2_1.fastq --pe2-2 '+'testdata_'+file2+'_bowtie2_2.fastq --pe3-1 '+'testdata_'+file3+'_bowtie2_1.fastq --pe3-2 '+'testdata_'+file3+'_bowtie2_2.fastq --pe4-1 '+'testdata_'+file4+'_bowtie2_1.fastq --pe4-2 '+'testdata_'+file4+'_bowtie2_2.fastq -o Assembly/'
    #run spades command above 
    os.system(runSPAdes)
    #write spades command to the log file
    with open('PipelineProject.log', 'a') as log_file:
        log_file.write(str(runSPAdes))
#call the SPAdes function with the provided input sequences 
SPAdes(sequences)

#define a function numberContigs
def numberContigs():
    #open file and write the contigs
    file = open('Step4Contigs.txt','w')
    #initialize a count of the contigs
    intial = 0

    #open the SPAdes assembly and go over its records
    file_path = SeqIO.parse('./Assembly/contigs.fasta','fasta')
    for i in file_path:
        #get the length of the current contigs sequence
        variable = len(i.seq)

        #if the contig is longer than 1000 bp, add the contig count and write the contig to the output file
        if variable > 1000:
            intial = intial+1

            file.write('>'+str(i.id) +'\n' +str(i.seq) + '\n')

    #close the file
    file.close()
    #write the contigs above 1000 in the log file
    with open('PipelineProject.log', 'a') as log_file:
        log_file.write('There are '+str(intial)+' contigs > 1000 bp in the assembly.\n')

#call the numberContigs function
numberContigs()

#define the function countnumberContigs
def countnumberContigs():
    #open the file that contains the contigs
    contigsFile = open('Step4Contigs.txt','r')

    #parse the contigs as a fasta file
    file = SeqIO.parse('Step4Contigs.txt','fasta')

    #create an empty list to store the length of each of the contigs
    countContigs = []

    #for loop to go through each contig in the file
    for i in file:
        #get the length of the contig
        variable = len(i.seq)
        #append the length to the list 
        countContigs.append(int(variable))

    #Calculate the total length of all the contigs in the file
    final = sum(countContigs)
    #close the contigs file
    contigsFile.close()
    #print the output into the log file
    with open('PipelineProject.log', 'a') as log_file:
        log_file.write('\n'+'There are '+ str(final)+'bp in the assembly.')
    
#call the function to run
countnumberContigs()

#define a function to get the longest contig
def contigLongest():
    #open the input file containing contigs and store it in the variable 'fileInput' 
    fileInput = open('Step4Contigs.txt','r')

    #parse the input file as fasta file using SeqIO and store it in the file
    file = SeqIO.parse('Step4Contigs.txt','fasta')

    #initialize max_length and max_record to none 
    max_length = None
    max_record = None

    #for loop to go over the record in the parsed fasta file
    for i in file:
        #this will get the length of the current contig sequence 
        variable = len(i.seq)

        #this will check if max_length' is equal to none and 'max_length' to check if the varibale is greater than it
        if (max_length ==None):
            #if so, then update the max_length and max_record with the current contig's length and sequence record 
            max_length = variable
            max_record = i
        #this will check if the varibale is greater than max_length 
        if variable > max_length:
            #if so, update the max_length with the current contig length and sequence
            max_length = variable
            max_record = i

    #close the input file
    fileInput.close()
    #write the sequence record of the longest contig to a new fasta file named 'contigLongest.fasta'
    SeqIO.write(max_record, 'contigLongest.fasta','fasta')

#call the function to run
contigLongest()

