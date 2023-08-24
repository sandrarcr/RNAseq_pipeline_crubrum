#################### Steps for creating a de novo transcriptome assembly for Corallium rubrum under ocean warming conditions ##################

#RNAseq pipeline

#Short reads
#N = 36

#TREATMENTS:
#control and high temperature (18C vs 25C) in   T0 
#                                               T1
#                                               T2
#                                               C0
#                                               C1
#                                               C2

########################################### Raw data quality control #############################################

#.md5 check | inside the directory where the unprocessed data is, do:

for file in *.md5; do md5sum -c $file; done

#HMGWJDSX3_4_86UDI-idt-UMI_1.fastq.gz: OK
#HMGWJDSX3_4_86UDI-idt-UMI_2.fastq.gz: OK

#Merge & rename files with individual ID.

#The rename function in marbits is util-linux rename (man7.org/linux/man-pages/man1/rename.1.html), which doesn't support s/.../.../ shit that you always do. 
#So, you have to do:

#rename 'part_you_want_to_change' 'thing_you_want_to_change_it_for' *extension_of_file(.fq.gz)

#files look like this:
#HMGWJDSX3_2_74UDI-idt-UMI_1.fastq.gz, so do:

rename 's/HMGWJDSX3_//' *fastq.gz
rename 's/HGWGKDSX3_//' *fastq.gz

# do it with HGWGKDSX3 too
#This is to remove the ID "HMGWJDSX3_1_" that mess up the concatenation

rename 's/UMI_1/UMI_R1/' *fastq.gz #to add which read is it. Do the same for R2
rename 's/UMI_2/UMI_R2/'

#Then, add to the /bin/$PATH the loop_merge.sh script:
#execute script

./loop_merge.sh

#loop_merge
for name in ./*.fastq.gz; 
do
    rnum=${name##*_}
    rnum=${rnum%%.*}

    sample=${name#*_}
    sample=${sample%%_*}

    cat "$name" >>"${sample}_$rnum.fastq.gz"
done

#what it does:
#This would iterate over all compressed Fastq files in the current directory and extract the sample name into the shell variable sample. 
#For all the filenames shown in the question, this would be XXUDI-idt-UMI_R1.
#The rnum variable will hold the R# bit at the end of the filename.
#The sample name is extracted by taking the filename and first removing everything up to and including the first _ character, and then removing everything after and including the first _ character from that result. The value for the rnum variable is extracted in a similar manner.
#The file is then simply appended onto the end of the aggregated file using cat >>. The output filename will be constructed from the sample name, the R#, and the string .fastq.gz. For the shown files, this will be XXUDI-idt-UMI_R1.fastq.gz.

#check number of reads for R1 and R2 with:

zcat 86UDI-idt-UMI_R1.fastq.gz |  echo $((`wc -l`/4))
58593553

#QUALITY CONTROL OF THE DATA:

module avail #command to check programs intalled in marbits --> verify you have what you need

#FastQC
#Trimmomatic
#BUSCO
#Transdecoder
#Trinity
#RSEM
#Salmon

########################################## 1) FastQC #############################################

#define variables:

echo $SLURM_ARRAY_TASK_ID

inputFile=$(awk "NR==$SLURM_ARRAY_TASK_ID" fastqFileList.txt)
echo $inputFile

#Design script: nano name_of_script.slurm

nano fastqc.slurm

#what's inside:

#a submission to marbits could be: It starts from #!/bin/bash

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=FASTQC_crubrum_rnaseq
##SBATCH --time=0-3:00              # Runtime in days-hours:minutes
#SBATCH --mem=6G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=2            # default = 1, max = 48
#SBATCH --output=FastQC.log.%A_%a.out     # File to which standard out will be written | path for output files
#SBATCH --error=FastQC.%A_%a.err      # File to which standard err will be written | provide path
#SBATCH --array=1-72%1                  #se pone un rango con el n de elementos que se tienen
##SBATCH --mail-type=ALL            # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=sandraramirezcalero@gmail.com  # Email to which notifications will be sent

#load modules here:

module load fastqc/0.11.7 #check for possible versions with:

# Get input file name from list
inputFile=$(awk "NR==$SLURM_ARRAY_TASK_ID" fastqFileList.txt)

# do all the interesting things
fastqc --outdir fastqc_output ${inputFile}

############ Script finish here ############

#execute
sbatch fastqc.slurm

# for f in /directory/of/the/data/*.fastq.gz; do fastqc --outdir /output/directory/fastqc_output $f ; done #the one i used in HK

#################################### 2) Trimmomatic #######################################
#To remove adapters that come from the illumina and to remove bad sequences

#define variables:

echo $SLURM_ARRAY_TASK_ID_1

inputFile1=$(awk "NR==$SLURM_ARRAY_TASK_ID_1" Ind_names_crubrum_R1.txt)
echo $inputFile1

echo $SLURM_ARRAY_TASK_ID_2

inputFile2=$(awk "NR==$SLURM_ARRAY_TASK_ID_2" Ind_names_crubrum_R2.txt)
echo $inputFile2

sbatch <script_name>

#the script may be:
#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=Trim_crubrum_rnaseq
#SBATCH --mem=6G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=2            # default = 1, max = 48
#SBATCH --output=trimmomatic.log.%A_%a.out     # File to which standard out will be written | path for output files
#SBATCH --error=trimmomatic.%A_%a.err      # File to which standard err will be written | provide path
#SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load trimmomatic/0.38 #check versions with module avail

# Get input file name from list
inputFile1=$(awk "NR==$SLURM_ARRAY_TASK_ID_1" Ind_names_crubrum_R1.txt)
inputFile2=$(awk "NR==$SLURM_ARRAY_TASK_ID_2" Ind_names_crubrum_R2.txt)

#run trimmomatic here:

for f1 in *_R1.fastq.gz
do
    f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
     trimmomatic PE -threads 10 -phred33 -trimlog f1_trimmolog $f1 $f2 "$f1"_paired.fq.gz "$f1"_unpaired.fq.gz "$f2"_paired.fq.gz "$f2"_unpaired.fq.gz ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10 LEADING:4 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 ;
done

trimmomatic PE -threads 10 -phred33 -trimlog f1_trimmolog ${inputFile1} ${inputFile2} "${inputFile1}"_paired.fq.gz "${inputFile1}"_unpaired.fq.gz ${inputFile2}_paired.fq.gz ${inputFile2}_unpaired.fq.gz ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10 LEADING:4 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

######## This script never worked, we needed to find a solution, so I did this:

tmux new-session #open an interactive session within marbits 
srun --pty /bin/bash #This will allocate an interactive session in a compute node and will asign to you 2 CPUs to work

#from the direct command line (and the respective directory where I wanted to run things, so I did:

for f1 in *_R1.fastq.gz
do
    f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
     trimmomatic PE -threads 10 -phred33 -trimlog f1_trimmolog $f1 $f2 "$f1"_paired.fq.gz "$f1"_unpaired.fq.gz "$f2"_paired.fq.gz "$f2"_unpaired.fq.gz ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10 LEADING:4 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 ;
done

#obtained:

CAS1_C1_R1.fastq.gz_paired.fq.gz # outputfiles from R1 
CAS1_C1_R1.fastq.gz_unpaired.fq.gz 

CAS1_C1_R2.fastq.gz_paired.fq.gz # outputfiles from R2
CAS1_C1_R2.fastq.gz_unpaired.fq.gz

# 'paired' output where both reads survived the processing --> these are the ones we keep
# 'unpaired' output where a read survived, but the partner read did not

#You can close the window but check how and you can leave the job running
tmux a #you can come back to your interactive session where the job is running

Ctrl + B and X #This will close entirely the interactive session and stop showing you have an active job in sqa

###################################### 3) FastQC after trimming #############################################

#within trimmomatic directory (where the outputfiles are), I created:

mkdir fastqc_output
#and provided the outputfile new names list as:

fastqFileList_trimmo.txt #where the names were: 

CAS_T1_2_R1.fastq.gz_paired.fq.gz
CAS_T1_2_R2.fastq.gz_paired.fq.gz

#then execute:

sbatch ../scripts/fastqc_trimmo.slurm

#and inside this:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=FASTQC_Trim_crubrum
##SBATCH --time=0-3:00              # Runtime in days-hours:minutes
#SBATCH --mem=6G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=5            # default = 1, max = 48
#SBATCH --output=FastQC.log.%A_%a.out     # File to which standard out will be written | path for output files
#SBATCH --error=FastQC.%A_%a.err      # File to which standard err will be written | provide path
#SBATCH --array=1-72%5                  #se pone un rango con el n de elementos que se tienen
##SBATCH --mail-type=ALL            # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=sandraramirezcalero@gmail.com  # Email to which notifications will be sent

#load modules here:

module load fastqc/0.11.7 #check for possible versions with:

# Get input file name from list
inputFile=$(awk "NR==$SLURM_ARRAY_TASK_ID" fastqFileList_trimmo.txt)

# do all the interesting things
fastqc --outdir fastqc_output ${inputFile}

############### Script ends here ################

############################################# 4) Kraken #################################################
#Taxonomic sequence classifier that assigns taxonomic labels to short DNA reads. It does this by examining the -mers within a read and querying a database with those K-mers. We need to eliminate bacteria, fungi and virus that don't belong to anything that it's coral.

##### 4.1) Kraken installation locally in my Marbits user ######

#Download and copy the kraken2-master folder from https://github.com/DerrickWood/kraken2
#go to folder where it was download and do:

#code
#!/bin/sh

#Download source code for kraken2 (NOTE: version 2.0.7-beta doesn't work)
download kraken2 version 2.0.8-beta version from https://github.com/DerrickWood/kraken2

cd '/home/sramirezc/programs/kraken2-master'# Go to the directory where kraken2 was extracted.

##  To install kraken2
./install_kraken2.sh '/home/sramirezc/programs/install_kraken2' # successful if you see the message: "Kraken 2 installation complete"

## Add the following directories to PATH - 

nano ~/.bashrc

## Add these to ~/.bashrc
export PATH=$PATH:/mnt/lustre/bio/users/sramirezc/programs/kraken2/kraken2
export PATH=$PATH:/mnt/lustre/bio/users/sramirezc/programs/kraken2/kraken2-build
export PATH=$PATH:/mnt/lustre/bio/users/sramirezc/programs/kraken2/kraken2-inspect
export PATH=$PATH:/mnt/lustre/bio/users/sramirezc/programs/kraken2/kraken2lib.pm
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

######## 4.2) Build Kraken2 custom database ###################

############### Step 4.2.1: Download the ncbi taxonomy names and tree and accession number to taxon maps
#download available databases for kraken in https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads, with:

wget -m ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/. 

./kraken2-build --download-taxonomy --db ../kraken2_database

############### Step 4.2.2: Install all required reference libraries (I used bacteria, fungus and virus)
./kraken2-build --download-library bacteria --db ../kraken2_db --no-masking #replace bacteria with fungi and viral
./kraken2-build --download-library fungi --db ../kraken2_db --no-masking
./kraken2-build --download-library viral --db ../kraken2_db --no-masking

############### Step 4.3: Build the database ##################

nano kraken-build-lib

#the script may be:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=kraken-build
#SBATCH --mem=6G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=12            # default = 1, max = 48
##SBATCH --output=c.log.%A_%a.out     # File to which standard out will be written | path for output files
##SBATCH --error=trimmomatic.%A_%a.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#run:

./kraken2-build --threads 12 --build --db ../kraken2_db/

############### Script ends here ################

sbatch kraken-build-lib

## Output:
#Database construction complete. [Total: 7h3m27.058s]

############### Step 4.4 - Classify my reads ##############################

#script 
nano kraken #Always check the paths and the name of the files before running:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=kraken-class
#SBATCH --mem=150G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=24            # default = 1, max = 48
##SBATCH --output=c.log.%A_%a.out     # File to which standard out will be written | path for output files
##SBATCH --error=trimmomatic.%A_%a.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#do interesting things here:

for R1 in *_R1.fastq.gz_paired.fq.gz; do R2=${R1/_R1.fastq.gz_paired.fq.gz}"_R2.fastq.gz_paired.fq.gz";output=${R1/_R1.fastq.gz_paired.fq.gz}"output.txt";report=${R1/_R1.fastq.gz_paired.fq.gz}"report.txt";classified=${R1/_R1.fastq.gz_paired.fq.gz}"classified";unclassified=${R1/_R1.fastq.gz_paired.fq.gz}"unclassified"; /home/sramirezc/sramirezc/programs/kraken2-master/kraken2 --db /home/sramirezc/sramirezc/programs/kraken2-master/kraken2_db/ --threads 24 --gzip-compressed --confidence 0.3 --output /home/sramirezc/sramirezc/crubrum_rnaseq/kraken/output/$output --report /home/sramirezc/sramirezc/crubrum_rnaseq/kraken/report/$report --paired --use-names --report-zero-counts --classified-out /home/sramirezc/sramirezc/crubrum_rnaseq/kraken/$classified#.fq --unclassified-out /home/sramirezc/sramirezc/crubrum_rnaseq/kraken/$unclassified#.fq $R1 $R2; bgzip /home/sramirezc/sramirezc/crubrum_rnaseq/kraken/*unclassified_1.fq; bgzip /home/sramirezc/sramirezc/crubrum_rnaseq/kraken/*unclassified_2.fq; done 

# change the directory folder if neccesary

sbatch kraken

#output

################################## 5) Design transcriptome reference using Trinity ################################

#script:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=Trinity_crubrum_rnaseq
#SBATCH --mem=150G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=Trinity.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=Trinity.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load trinity/2.13.2

# do the interesting things here

Trinity --seqType fq --max_memory 20G --left CAS1_C1_unclassified_1.fq.gz,CAS_C0_1_unclassified_1.fq.gz,CAS_C2_1_unclassified_1.fq.gz,CAS_T0_1_unclassified_1.fq.gz,CAS_T1_2_unclassified_1.fq.gz,CAS_T2_1_unclassified_1.fq.gz,LOP2_C1_unclassified_1.fq.gz,LOP_C0_1_unclassified_1.fq.gz,LOP_C2_1_unclassified_1.fq.gz,LOP_T0_1_unclassified_1.fq.gz,LOP_T1_1_unclassified_1.fq.gz,LOP_T2_1_unclassified_1.fq.gz --right CAS1_C1_unclassified_2.fq.gz,CAS_C0_1_unclassified_2.fq.gz,CAS_C2_1_unclassified_2.fq.gz,CAS_T0_1_unclassified_2.fq.gz,CAS_T1_2_unclassified_2.fq.gz,CAS_T2_1_unclassified_2.fq.gz,LOP2_C1_unclassified_2.fq.gz,LOP_C0_1_unclassified_2.fq.gz,LOP_C2_1_unclassified_2.fq.gz,LOP_T0_1_unclassified_2.fq.gz,LOP_T1_1_unclassified_2.fq.gz,LOP_T2_1_unclassified_2.fq.gz --CPU 30 --output /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/trinity_denovo_c.rubrum

### Script ends here ###

sbatch Trinity.slurm

#Check for Trinity.err inside //home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/

#output: trinity_denovo_c.rubrum.Trinity.fasta
#TOTAL TRANSCRIPTS FROM RAW DE NOVO TRANSCRIPTOME: 533146

#Now, we need to check we didn't end up with a transcriptome that actually has the genes that come from the Coral, so we do:

########################### 6) Transcriptome Assembly Quality Assessment #########################
#mapping/align with Bowtie2
#see completeness with BUSCO
#use N50
#transdecoder

################### 6.1) Transdecoder ########################
#identifies candidate coding regions within transcript sequences because it is important to know that the contigs created during de novo assembly actually transcribed real genes and to avoid reduncy as you end up with A LOT of data.

#script:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_transdecoder
#SBATCH --mem=10G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=transdecoder.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=transdecoder.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load transdecoder/5.5.0
module load perl/5.28

# do the interesting things here

TransDecoder.LongOrfs -t crubrum.trinity.fasta --output_dir /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/transdecoder/

#### Script ends here :) #####

#We obtain longest_orfs.pep: 204254 longest ORFs

########## Step 6.1.1) Including homology searches as ORF retention criteria ##########

#script

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_blastp_uniprot
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=blastp_uniprot.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastp_uniprot.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1

# do the interesting things here:

blastp -query longest_orfs.pep -db /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/databases/uniprot/uniprot_021222.fasta -out ./c.rubrum_blastp_uniprot.outfmt6.txt -evalue 1e-5 -num_threads 30 -max_target_seqs 1 -outfmt 6

#################### Script ends here ###################

#longest_orfs.pep   : all ORFs meeting the minimum length criteria, regardless of coding potential.
#longest_orfs.gff3  : positions of all ORFs as found in the target transcripts
#longest_orfs.cds   : the nucleotide coding sequence for all detected ORFs
#transcripts.fasta.transdecoder.pep : peptide sequences for the final candidate ORFs; all shorter candidates within longer ORFs were removed.
#transcripts.fasta.transdecoder.cds  : nucleotide sequences for coding regions of the final candidate ORFs
#transcripts.fasta.transdecoder.gff3 : positions within the target transcripts of the final selected ORFs
#transcripts.fasta.transdecoder.bed  : bed-formatted file describing ORF positions, best for viewing using GenomeView or IGV.

############ Step 6.1.2) Predict the likely coding regions ##########

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_transdecoder
#SBATCH --mem=10G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=20            # default = 1, max = 48
#SBATCH --output=transdecoder_pre.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=transdecoder_pre.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load transdecoder/5.5.0
module load perl/5.28

# do the interesting things here

TransDecoder.Predict -t ../../crubrum.trinity.fasta --retain_blastp_hits c.rubrum_blastp_uniprot.outfmt6.txt --single_best_only --output_dir /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/transdecoder/

#### Script ends here :) #####

#You obtain several files, but focus on: crubrum.trinity.fasta.transdecoder.pep and crubrum.trinity.fasta.transdecoder.cds and it looks like this:

>TRINITY_DN0_c0_g1_i2.p1 TRINITY_DN0_c0_g1~~TRINITY_DN0_c0_g1_i2.p1  ORF type:complete len:857 (-),score=161.30,sp|Q6ZMV9|KIF6_HUMAN|49.424|0.0 TRINITY_DN0_c0_g1_i2:1108-3678(-)
MVKYGIQIFARVKPTRGKTGDYDCEDEDDGFSYVSFNVPKDLARDFINNKKEIYKFRFNK
AFEKDIKQDDVFQYVAKGVIDNVLSGYNGTIFAYGQTGSGKTFTITGGAERYADRGIIPR
TLSYMFECFEKNPESVYTSHVSYLEIYNENGYDLLDPKHEACKLEDLPKIALMEDNDGNI
HLKNLSLHQANNEEEALNWLFLGDTNRMIAETPMNQASTRSHCIFTMHVSSREPGSATLR
RAKLHLVDLAGSERIHKSNIDGTLLTEAKYINLSLHYLEQVIVALSEKSRTHIPYRNSLL
TSVLRDSLGGNCRTTMIATLSIDRKNIDESISTCRFAQRVALIKNDAILNEEIDPKLMII
KLKQEIQQLKDELAISSGQEYRGELTDEDIERLKFMIKAYIEDRDPEALLSVGADMRKIN
LSFKLFKGHVLERKANSAPAAIRSSLTSQDSSYLSDQNESKKLKELIQQRDNEINILVGM
LKRERERVAAGLKTVDDRTWDHRDERRVPKKGVSRQSSSSSTHGESSEHQSQTTLEERGQ
PQAIKEEAGEGHNRKKISHSKQLADDMSQGRREAFEVFRSDYPHNATIEENKRTLKQRYS
EAKALGEQVNNSRKRINTIKTQIEQHRIRRSMHVNGDEDNIDEEDRLRNSIEEEKSRYKD

#Total transcripts with longest ORF and homology: 124373

############# Step 6.1.3) Get complete transcript seq for final selected contigs. It selects the column that says gene

awk 'BEGIN{FS="\t"}{if ($3 == "gene") print}' *.transdecoder.gff3 > final_assembly_after_transdecoder_genes.gff3 

cut -f 1,4,5 *transdecoder_genes.gff3 > crubrum_final_transdecoder_genes.bed #cut columns 1, 4 and 5 to the output

bedtools getfasta -fi crubrum.trinity.fasta -bed crubrum_final_transdecoder_genes.bed -fo crubrum_final_reference.fasta #using bedtools to extract the nucleotide sequences corresponding to the genomic coordinates specified in the crubrum_final_transdecoder_genes.bed file, and write them to a new file called crubrum_final_reference.fasta!!! Keep it for further analysis.

less crubrum_final_reference.fasta
# This is the definitive

######################### 6.2) Bowtie2 ####################################

######################### 6.2.1) Bowtie2_index_build ##############################
#Mapping against your reads

#scrcipt
#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=c.rubrum_bowtie2_index
#SBATCH --mem=6G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=Bowtie2-index.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=Bowtie2-index.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load bowtie2/2.3.4.1

# do the interesting things here

bowtie2-build crubrum_final_reference.fasta crubrum_bowtie_index_2

### Script ends here #####

sbatch bowtie2_index.slurm

######################### 6.2.2 Bowtie2 alignment ##########################

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=c.rubrum_bowtie2_align
#SBATCH --mem=10G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=Bowtie2-align.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=Bowtie2-align.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load bowtie2/2.3.4.1

# do the interesting things here

for R1 in *_R1.fq.gz; do R2=${R1/R1.fq.gz}"R2.fq.gz"; output1=${R1/_R1.fq.gz}"_stats"; bam=${R1/_R1.fq.gz}".bam"; bowtie2 -q -p 30 -x crubrum_bowtie_index --met-file /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/Bowtie2/$output1 -1 $R1 -2 $R2 2>align_stats_$output1 | samtools view -@10 -Sb -o /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/Bowtie2/Bam/$bam ; done

#### Script ends here :) #####

#Align_stats of one individual

26451563 reads; of these:
  26451563 (100.00%) were paired; of these:
    9515936 (35.97%) aligned concordantly 0 times
    4224328 (15.97%) aligned concordantly exactly 1 time
    12711299 (48.06%) aligned concordantly >1 times
    ----
    9515936 pairs aligned concordantly 0 times; of these:
      641133 (6.74%) aligned discordantly 1 time
    ----
    8874803 pairs aligned 0 times concordantly or discordantly; of these:
      17749606 mates make up the pairs; of these:
        12486770 (70.35%) aligned 0 times
        494202 (2.78%) aligned exactly 1 time
        4768634 (26.87%) aligned >1 times
76.40% overall alignment rate

#"Alignment" is the process by which we discover how and where the read sequences are similar to the reference sequence. An "alignment" is a result from this process, specifically: an alignment is a way of "lining up" some or all of the characters in the read with some characters from the reference in a way that reveals how they're similar.

  Read:      GACTGGGCGATCTCGACTTCG
             |||||  |||||||||| |||
  Reference: GACTG--CGATCTCGACATCG


#We use alignment to make an educated guess as to where a read originated with respect to the reference genome/transcriptome. 

#By default, Bowtie 2 performs end-to-end read alignment. That is, it searches for alignments involving all of the read characters.

Read:      GACTGGGCGATCTCGACTTCG
Reference: GACTGCGATCTCGACATCG

Alignment:
 
  Read:      GACTGGGCGATCTCGACTTCG
             |||||  |||||||||| |||
  Reference: GACTG--CGATCTCGACATCG

#An alignment score quantifies how similar the read sequence is to the reference sequence aligned to. The higher the score, the more similar they are. A score is calculated by subtracting penalties for each difference (mismatch, gap, etc.)
#The best possible alignment score in end-to-end mode is 0, which happens when there are no differences between the read and the reference.

#A "paired-end" or "mate-pair" read consists of pair of mates, called mate 1 and mate 2. Pairs come with a prior expectation about (a) the relative orientation of the mates, and (b) the distance separating them on the original DNA molecule. Exactly what expectations hold for a given dataset depends on the lab procedures used to generate the data. 

#A pair that aligns with the expected relative mate orientation and with the expected range of distances between mates is said to align "concordantly". 
#If both mates have unique alignments, but the alignments do not match paired-end expectations (i.e. the mates aren't in the expected relative orientation, or aren't within the expected distance range, or both), the pair is said to align "discordantly".

Mate 1:    GCAGATTATATGAGTCAGCTACGATATTGTT
Mate 2:                               TGTTTGGGGTGACACATTACGCGTCTTTGAC
Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

############### 6.3) Trinity Transcript Quantification #############

#remove from headers in crubrum_final_reference.fasta:

sed -i 's/:.*$//' crubrum_final_reference.fasta

#from this: >TRINITY_DN0_c0_g1_i2:1-3918 
#to this: >TRINITY_DN0_c0_g1_i2

######### 6.3.1) create reference:

#script:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_RSEM_ref
#SBATCH --mem=10G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=RSEM_ref.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=RSEM_ref.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load trinity/2.13.2
module load bowtie2/2.3.4.1
module load salmon/1.8.0
module load rsem/1.3.3
module load perl/5.28

# do the interesting things here

align_and_estimate_abundance.pl  --transcripts crubrum_final_reference.fasta --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference

#### Script ends here :) #####

less RSEM_ref.err:

#CMD: touch /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta.bowtie2.started
#CMD: bowtie2-build /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta.bowtie2
#Building a SMALL index
#WARNING - appears that another process has started the rsem-prep step... proceeding with caution.
#CMD: touch /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta.RSEM.rsem.prepped.started
#CMD: rsem-prepare-reference  --transcript-to-gene-map /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta.gene_trans_map /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/RSEM_out/RSEM_ref_est/crubrum_final_reference.fasta.RSEM
#Only prepping reference. Stopping now.
######## Script ends here ######

###### 6.3.2) using align_and_estimate_abundance.pl ########

#RSEM
#Bowtie2

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_RSEM
#SBATCH --mem=150G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=48            # default = 1, max = 48
#SBATCH --output=RSEM_est.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=RSEM_est.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load trinity/2.13.2
module load bowtie2/2.3.4.1
module load rsem/1.3.3
module load perl/5.28

# do the interesting things here

for R1 in *_R1.fq.gz; do R2=${R1/R1.fq.gz}"R2.fq.gz"; output=${R1/_R1.fq.gz}"_RSEM_out"; align_and_estimate_abundance.pl  --transcripts crubrum_final_reference.fasta --seqType fq --left $R1 --right $R2 --est_method RSEM --aln_method bowtie2 --output_dir $output --SS_lib_type RF --gene_trans_map crubrum_final_reference.fasta.gene_trans_map --thread_count 48; done

##### Script ends here ####

###################### 6.4) Build Transcript and Gene Expression Matrices ##########

#### You get these files per sample (R1&R2):
#Folders: _RSEM_out 
#files: RSEM.stat, bowtie2.bam bowtie2.bam.for_rsem.bam, bowtie2.bam.ok, RSEM.isoforms.results.ok, RSEM.genes.results, RSEM.isoforms.results

#The most important are RSEM.genes.results and RSEM.isoforms.results, but you have to change the name of each to build the matrix as follows. Use the isoforms. results file always

for directory in *unclassified_RSEM_out*; do
  pushd "$directory"
  
  for filename in *results; do
    extension="${filename}"
    
    
    target_filename="${directory}_${extension}"
    
    mv "$filename" "${target_filename}"
    
  done
  popd
done

#and then put all isoforms.results files in one folder and change what is left like this:
#Do this in the folder where the *_RSEM_out folders of all samples are

mkdir isoforms.results #then, do

cp ./*/*isoforms.results ./isoforms.results/ #this says from the present directory where there are more directories with files named *isoforms.results, copy them to ./isoforms.results

#result:

cd isoforms.results

CAS1_C1_unclassifiedRSEM_out_RSEM.isoforms.results   CAS_T1_2_unclassifiedRSEM_out_RSEM.isoforms.results  LOP_C2_1_unclassifiedRSEM_out_RSEM.isoforms.results
CAS2_C1_unclassifiedRSEM_out_RSEM.isoforms.results   CAS_T1_3_unclassifiedRSEM_out_RSEM.isoforms.results  LOP_C2_2_unclassifiedRSEM_out_RSEM.isoforms.results
CAS3_C1_unclassifiedRSEM_out_RSEM.isoforms.results   CAS_T1_4_unclassifiedRSEM_out_RSEM.isoforms.results  LOP_C2_3_unclassifiedRSEM_out_RSEM.isoforms.results
CAS_C0_1_unclassifiedRSEM_out_RSEM.isoforms.results  CAS_T2_1_unclassifiedRSEM_out_RSEM.isoforms.results  LOP_T0_1_unclassifiedRSEM_out_RSEM.isoforms.results

#rename

rename '_unclassifiedRSEM_out_RSEM' '' *isoforms.results

CAS1_C1.isoforms.results   CAS_C2_2.isoforms.results  CAS_T1_4.isoforms.results             LOPB_CA.isoforms.results   LOP_T0_1.isoforms.results  LOP_T2_2.isoforms.results
CAS2_C1.isoforms.results   CAS_C2_3.isoforms.results  CAS_T2_1.isoforms.results             LOP_C0_1.isoforms.results  LOP_T0_2.isoforms.results  LOP_T2_3.isoforms.results
CAS3_C1.isoforms.results   CAS_T0_1.isoforms.results  CAS_T2_2.isoforms.results             LOP_C0_2.isoforms.results  LOP_T0_3.isoforms.results

#then run the following:

###################### 6.5) Build Transcript and Gene Expression Matrices ###################

#script

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_RSEM_matrix
#SBATCH --mem=10G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=10            # default = 1, max = 48
#SBATCH --output=RSEM_matrix.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=RSEM_matrix.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load trinity/2.13.2
module load bowtie2/2.3.4.1
module load rsem/1.3.3
module load perl/5.28

# do the interesting things here

abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/downstream_analyses/RSEM/crubrum_final_reference.fasta.gene_trans_map --out_prefix crubrum CAS1_C1.isoforms.results CAS2_C1.isoforms.results CAS3_C1.isoforms.results CAS_C0_1.isoforms.results CAS_C0_2.isoforms.results CAS_C0_3.isoforms.results CAS_C2_1.isoforms.results CAS_C2_2.isoforms.results CAS_C2_3.isoforms.results CAS_T0_1.isoforms.results CAS_T0_2.isoforms.results CAS_T0_3.isoforms.results CAS_T1_2.isoforms.results CAS_T1_3.isoforms.results CAS_T1_4.isoforms.results CAS_T2_1.isoforms.results CAS_T2_2.isoforms.results CAS_T2_3.isoforms.results LOP2_C1.isoforms.results LOP3_C1.isoforms.results LOPB_CA.isoforms.results LOP_C0_1.isoforms.results LOP_C0_2.isoforms.results LOP_C0_5.isoforms.results LOP_C2_1.isoforms.results LOP_C2_2.isoforms.results LOP_C2_3.isoforms.results LOP_T0_1.isoforms.results LOP_T0_2.isoforms.results LOP_T0_3.isoforms.results LOP_T1_1.isoforms.results LOP_T1_2.isoforms.results LOP_T1_3.isoforms.results LOP_T2_1.isoforms.results LOP_T2_2.isoforms.results LOP_T2_3.isoforms.results

##### Script ends here ####

#Output files

#-rw-r--r--  1 sramirezc bio  60934425 Dec  5 10:00 crubrum.gene.counts.matrix
#-rw-r--r--  1 sramirezc bio  70160052 Dec  5 10:02 crubrum.gene.TMM.EXPR.matrix
#-rw-r--r--  1 sramirezc bio  34563094 Dec  5 10:00 crubrum.gene.TPM.not_cross_norm
#-rw-r--r--  1 sramirezc bio       524 Dec  5 10:01 crubrum.gene.TPM.not_cross_norm.runTMM.R
#-rw-r--r--  1 sramirezc bio      1840 Dec  5 10:02 crubrum.gene.TPM.not_cross_norm.TMM_info.txt
#-rw-r--r--  1 sramirezc bio 111106892 Dec  5 10:00 crubrum.isoform.counts.matrix
#-rw-r--r--  1 sramirezc bio 128812305 Dec  5 10:01 crubrum.isoform.TMM.EXPR.matrix
#-rw-r--r--  1 sramirezc bio  63721483 Dec  5 10:00 crubrum.isoform.TPM.not_cross_norm
#-rw-r--r--  1 sramirezc bio       530 Dec  5 10:00 crubrum.isoform.TPM.not_cross_norm.runTMM.R
#-rw-r--r--  1 sramirezc bio      1846 Dec  5 10:01 crubrum.isoform.TPM.not_cross_norm.TMM_info.txt

#but the most important are:

crubrum.gene.counts.matrix --> the estimated RNA-seq fragment counts (raw counts)
crubrum.isoform.TPM.not_cross_norm --> a matrix of TPM expression values (not cross-sample normalized)
crubrum.isoform.TMM.EXPR.matrix --> a matrix of TMM-normalized expression values 

#he 'counts.matrix' file is used for downstream analyses of differential expression. The TMM.EXPR.matrix file is used as the gene expression matrix in most other analyses. For information on the importance of TMM (or cross-sample normalization in general)

#After this, you have to normalize your data and filter lowly expressed transcripts. Do this:

#script
#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_fil_matrix
#SBATCH --mem=10G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=10            # default = 1, max = 48
#SBATCH --output=RSEM_fil.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=RSEM_fil.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

#load modules here:

module load trinity/2.13.2
module load bowtie2/2.3.4.1
module load rsem/1.3.3
module load perl/5.28

# do the interesting things here

filter_low_expr_transcripts.pl --matrix crubrum.gene.TPM.not_cross_norm --transcripts /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/crubrum.trinity.fasta --min_expr_any 10 --highest_iso_only --trinity_mode

##### Script ends here ####

#Matrix is filtered.

##################### 6.6) Counting Full Length Trinity Transcripts with BLAST ##################

#One metric for evaluating the quality of a transcriptome assembly is to examine the number of transcripts that were assembled that appear to be full-length or nearly full-length.
# For non-model organisms, no such reference transcript set is available. If a high quality annotation exists for a closely related organism, then one might compare the assembled transcripts to that closely related transcriptome to examine full-length coverage. In other cases, a more general analysis to perform is to align the assembled transcripts against all known proteins and to determine the number of unique top matching proteins that align across more than X% of its length.

#Dowload from SwissProt (ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) in https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts
#use filezilla as the application to download uniprot database

#copy to marbits
#open a tunnel, then do:

 scp -r uniprot/ sramirezc@161.111.137.246:/home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/databases/

#or

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

#We have to chop the header cuz it'll mess 

zcat uniprotkb_swissprot.gz | awk '{if (/^>/) { print ">" $2} else { print $_}}' > swissprot.fa

#With this, we pass from this: 

zless uniprot_sprot.fasta.gz

>sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1
MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS
EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD
AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL
EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD
SFRKIYTDLGWKFTPL

#To this:le

less swissprot.fa

>Q6GZX4
MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPS
EKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLD
AKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHL
EKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDD
SFRKIYTDLGWKFTPL

####################################### blastx ############################

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_blastx
#SBATCH --mem=20G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=blastx_swissprot.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_swissprot.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1

# do the interesting things here:

blastx -query /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/trinity_denovo_c.rubrum.Trinity.fasta -db /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/databases/swissprotdb/swissprot.fa -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/swissprot_res/trinity_denovo_c.rubrum_blastx.outfmt6.txt -max_hsps 1 -evalue 1e-5 -num_threads 30 -max_target_seqs 5 -outfmt 6

############################## Script ends here ###########################

#output

#TRINITY_DN109637_c0_g1_i1       Q09575  29.386  228     150     9       669     4       307     529     1.22e-08        58.5
#TRINITY_DN109600_c0_g1_i1       P08764  54.545  88      39      1       354     94      455     542     3.25e-23        96.7
#TRINITY_DN109600_c0_g1_i1       P25241  50.526  95      40      2       345     82      455     549     5.25e-20        87.4
#TRINITY_DN109600_c0_g1_i1       P40815  49.474  95      41      2       345     82      453     547     1.91e-19        85.9
#TRINITY_DN109603_c0_g1_i1       Q9FD71  37.327  217     134     1       756     112     74      290     8.61e-49        168

#examine the percent of the target being aligned to by the best matching Trinity transcript, like so:

#script

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_blastx
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=blastx_swissprot_add.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_swissprot_add.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1
module load trinity/2.13.2

# do the interesting things here:

analyze_blastPlus_topHit_coverage.pl trinity_denovo_c.rubrum_blastx.outfmt6.txt ../../../c.rubrum_reference/crubrum.trinity.fasta ../databases/swissprotdb/swissprot.fa

############## Script ends here

#output: trinity_denovo_c.rubrum_blastx.outfmt6.txt.w_pct_hit_length, so do:

less trinity_denovo_c.rubrum_blastx.outfmt6.txt.w_pct_hit_length

#qseqid sseqid  pident  length  mismatch        gapopen qstart  qend    sstart  send    evalue  bitscore        db_hit_len      pct_hit_len_aligned     hit_descr
TRINITY_DN28918_c0_g1_i1        Q7Z4Q2  33.459  266     170     4       33      827     1       260     3.86e-26        117     260     38.24
TRINITY_DN48562_c0_g1_i10       P05099  35.882  170     95      5       727     1218    271     432     3.20e-13        76.3    162     32.86
TRINITY_DN3581_c0_g1_i5 B2RYN7  47.509  562     254     5       1010    2671    49      577     8.39e-153       471     529     91.05
TRINITY_DN19560_c0_g1_i1        P20911  78.369  282     61      0       38      883     14      295     8.63e-153       443     282     80.11
TRINITY_DN30833_c0_g1_i2        O15254  52.616  688     319     6       134     2185    16      700     0.0     739     685     97.86

#this adds fields to the blastx output file, including the top hit's length, percent of the hit's length included in the alignment to the Trinity transcript, and the header description for that database entry.

#hit_pct_cov_bin        count_in_bin    >bin_below
100     5651    5651
90      2376    8027
80      1753    9780
70      1668    11448
60      1840    13288
50      2148    15436
40      2845    18281
30      3762    22043
20      5099    27142
10      2963    30105

#There are 2376 proteins that each match a Trinity transcript by >80% and <= 90% of their protein lengths.
#There are 8027 proteins that are represented by nearly full-length transcripts, having >80% alignment coverage.
#There are 5651 proteins that are covered by more than 90% of their protein lengths.

############# #To obtain a rough list of genes to see what it's inside de assembly, i did this:

#I did another blastx using the uniprot database downloaded from https://www.uniprot.org/uniprotkb?query=%2A with only reviewed records (568363)
#Built the database as I did before but did not chop the ID. Then I did this:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_blastx_uniprot
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=blastx_uniprot.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_uniprot.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen


###### load modules

module load blast/2.7.1

# do the interesting things here:

blastx -query /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/crubrum.trinity.fasta -db /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/databases/uniprot/uniprot_021222.fasta -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/uniprot_res/c.rubrum_blastx_uniprot.outfmt6.txt -max_hsps 1 -evalue 1e-5 -num_threads 30 -max_target_seqs 5 -outfmt 6

#################### Script ends here ###################3

#OUTPUT FILE:

c.rubrum_blastx_uniprot.outfmt6.txt

##I input the following command to get the % of the target being alinged to the best matching at trinity transcripts:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_blastx
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=blastx_swissprot_add.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_swissprot_add.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1
module load trinity/2.13.2

# do the interesting things here:

analyze_blastPlus_topHit_coverage.pl c.rubrum_blastx_uniprot.outfmt6.txt ../../../c.rubrum_reference/crubrum.trinity.fasta  ../databases/uniprot/uniprot_021222.fasta

### Script ends here

####

#hit_pct_cov_bin        count_in_bin    >bin_below
100     5665    5665
90      2378    8043
80      1764    9807
70      1676    11483
60      1847    13330
50      2160    15490
40      2845    18335
30      3791    22126
20      5120    27246
10      3009    30255


Table 1 Number of ESTs used for comparative analysis of each species and number of BLASTP hits between the red coral transcrip- tome and the EST bank of each species (e value = 10e10)

######################## 6.7) Annotation ######################

#Here we also include what we just did: Swissprot/Uniprot

#We want to know which transcripts have an accurate gene annotation to predict their functions. Since octocorals are cnidarians poorly studied, we have to be creative and use a highly known curated database such as swissprot (the one we did before), and create a couple of databases more with available fasta annotations from other cnidarians and corals. So we do this:

#Went to NCBI and downloaded a fasta containing the genes from softcorals only:

https://www.ncbi.nlm.nih.gov/nuccore:

#Corallium rubrum (658)
#Hemicorallium imperiale (223)
#Hemicorallium laauense (194)
#Corallium japonicum (120)
#Pleurocorallium niveum (102)
#Pleurocorallium secundum (82)
#Pleurocorallium konojoi (80)
#Pleurocorallium elatius (64)
#Corallium tortuosum (56)
#Pleurocorallium porcellanum (55)
#Hemicorallium ducale (52)
#Corallium sp. 6 THT-2013 (39)
#Pleurocorallium thrinax (36)
#Pleurocorallium borneense (34)
#Hemicorallium regale (29)
#Paracorallium sp. 4 THT-2013 (27)
#Hemicorallium niobe (23)
#Hemicorallium abyssale (22)
#Pleurocorallium inutile (21)
#Corallium cf. elatius NEA-2012 (20)

#Total of 2303 record looking like this:

>MT558779.1 Coralliidae sp. HI-008 16S ribosomal RNA gene, partial sequence; and NADH dehydrogenase subunit 2 gene, partial cds; mitochondrial
CAGCTCCGTTTCTATCTACAATGTTCAATCCCAACTTTTCTGAGTTCTAGTACGCAAGGAATGGTCTCAG
CGATGCTCGACTATATTGGCCCCATCGACGTAATCAACTTCTGGCTGCTGCAAAGAAGGAGAACAAAAGG
TTTGGGATTAATAGGGTGCCACAAAAGTGGCGCACCTTTAGGTGCCATGTGGGTGCATAGCCCCTGGCAC
CCTATGGAATTAGCACTAGGGCTTATTATACTAATAGTGCTAACGTATGGATTGAAGGCCCCGGCTCTAA
GATTAGCAATATTCTGTTTGGGAGTGACACTACTTATGGGTGCTGCTGGGTTATTAGCTGAGCCCCATCT
GCTATGTTGGACACAGGCTATTAAGATGTTGGTGATGCTAAGTGGGTTAGCC

############ 6.7.1) make database with downloaded fasta file --custom databases ##############

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=softcoral_blastdb
#SBATCH --mem=100G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=48            # default = 1, max = 48
#SBATCH --output=blastdb_softcoral.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastdb_softcoral.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1

# do the interesting things here:

makeblastdb -in softcoralgenes.fasta -dbtype nucl -parse_seqids
#nucl because we're gonna use it with crubrum_reference

###blastx -query /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/crubrum.trinity.fasta -db /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/databases/uniprot/uniprot_021222.fasta -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/uniprot_res/c.rubrum_blastx_uniprot.outfmt6.txt -max_hsps 1 -evalue 1e-5 -num_threads 30 -max_target_seqs 5 -outfmt 6

#run blastn since it's nucleotide:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum-soft_blastx
#SBATCH --mem=100G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=48            # default = 1, max = 48
#SBATCH --output=blastx_soft.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_soft.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1

# do the interesting things here:

blastn -query /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/crubrum.trinity.fasta -db /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/databases/softcorals/softcoralgenes.fasta -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/softcoral_res/c.rubrum_blastx_softcorals.outfmt6.txt -max_hsps 1 -evalue 1e-5 -num_threads 48 -max_target_seqs 5 -outfmt 6

#################### Script ends here ###################

################### 6.7.2) Other corals: #####################

#Utilizando genes anotados de corales de Ensembl:
#-rw-r--r-- 1 sramirezc bio 76M Dec 29 12:25 Acropora_millepora-GCA_013753865.1-2022_03-cds.fa
#-rw-r--r-- 1 sramirezc bio 58M Dec 29 12:25 Dendronephthya_gigantea-GCA_004324835.1-2021_11-cds.fa
#-rw-r--r-- 1 sramirezc bio 55M Dec 29 12:25 Orbicella_faveolata-GCA_002042975.1-2021_12-cds.fa
#-rw-r--r-- 1 sramirezc bio 45M Dec 29 12:25 Pocillopora_damicornis-GCA_003704095.1-2021_11-cds.fa
#-rw-r--r-- 1 sramirezc bio 60M Dec 29 12:25 Stylophora_pistillata-GCA_002571385.1-2021_11-cds.fa

#concatenate all fasta downloaded: https://rapid.ensembl.org/info/about/species.html, filtering by "coral"

cat Acropora_millepora-GCA_013753865.1-2022_03-cds.fa Dendronephthya_gigantea-GCA_004324835.1-2021_11-cds.fa Orbicella_faveolata-GCA_002042975.1-2021_12-cds.fa Pocillopora_damicornis-GCA_003704095.1-2021_11-cds.fa Stylophora_pistillata-GCA_002571385.1-2021_11-cds.fa > corals_emble.fasta

#modify ID in order to creat blast db:

#current ID:
>XM_029343860.2 cds Amil_v2.1:CM024249.1:23834324:23853178:-1 gene:LOC114964514 gene_biotype:protein_coding transcript_biotype:protein_coding description:"multiple C2 and transmembrane domain-containing protein 1-like"

#we conserve the >XM...(accession n) etc at the beginning, but we also want the name of the gene. So we do:
perl -pe 's/ cds (.+?)description:(.*)/_$2/' corals_emble.fa > corals_chopped.fa
#output:

>XM_029339785.2_"SNARE-associated protein Snapin-like" #blast doesn't like spaces, so we do:

sed 's|\ |_|g' corals_chopped.fa > corals_chopped1.fa #to replace every ' ' (space) with '_' (underscore)

#output:

 head corals_chopped.fa
>XM_029339785.2_"SNARE-associated_protein_Snapin-like"
ATGGCCGCCGAGGTAAAGGTAAAGAAAGATGCATTTGAGCAAGGAAATCCGACCGCTGTA
GAAAATCAGAACGCACTAGCAGAAGGAATCATGCAGATATTCAAACCAGCAGTCGAGGAA
CTTGATGACAAAGTCTTAAACGTCAGACAAAGCCAAGTTGAACTTAGAGAACAAATTGAC
AAGCTTTCTCAAGATCTTCACAGGTTATCTGAACTTCAAGAAATTCCAGTTGACCTGGAG

#run makeblastdb

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=coral_blastdb
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=40            # default = 1, max = 48
#SBATCH --output=blastdb_coralembl.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastdb_coralembl.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1

# do the interesting things here:

makeblastdb -in corals_chopped1.fa -dbtype nucl -parse_seqids -out coralsdb

#nucl because we're gonna use it with crubrum_reference

#################### Script ends here ###################

#output:

-rw-r--r-- 1 sramirezc bio  16M Dec 30 14:07 coralsdb.nhr
-rw-r--r-- 1 sramirezc bio 1.9M Dec 30 14:07 coralsdb.nin
-rw-r--r-- 1 sramirezc bio 632K Dec 30 14:07 coralsdb.nog
-rw-r--r-- 1 sramirezc bio  22M Dec 30 14:07 coralsdb.nsd
-rw-r--r-- 1 sramirezc bio 382K Dec 30 14:07 coralsdb.nsi
-rw-r--r-- 1 sramirezc bio  64M Dec 30 14:07 coralsdb.nsq

#They're nucleotides

#Now, run blastn

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=coralsdb_blastn
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=blastn_corals.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastn_corals.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

###### load modules

module load blast/2.7.1

# do the interesting things here:

blastn -query /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/crubrum.trinity.fasta -db coralsdb -out c.rubrum_corals.outfmt6.txt -max_hsps 1 -evalue 1e-5 -num_threads 30 -max_target_seqs 5 -outfmt 6

#################### Script ends here ###################

#output

head c.rubrum_corals.outfmt6.txt

TRINITY_DN108947_c0_g1_i1 XM_028544344.1_"uncharacterized_LOC114523428,_transcript_variant_X177.714   175     31      6       3       170     1039    1212    2.75e-20        100
TRINITY_DN108947_c0_g1_i1 XM_028544339.1_"uncharacterized_LOC114523428,_transcript_variant_X177.714   175     31      6       3       170     1039    1212    2.75e-20        100
TRINITY_DN108947_c0_g1_i1 XM_028544345.1_"uncharacterized_LOC114523428,_transcript_variant_X177.714   175     31      6       3       170     1039    1212    2.75e-20        100
TRINITY_DN108947_c0_g1_i1 XM_028544341.1_"uncharacterized_LOC114523428,_transcript_variant_X177.714   175     31      6       3       170     1039    1212    2.75e-20        100

################# 6.5) N50 ######################

module load trinity/2.13.2 #then, do:

TrinityStats.pl trinity_denovo_c.rubrum.Trinity.fasta > trinity_N50.txt

#output

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  293131
Total trinity transcripts:      533146
Percent GC: 40.05

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3669
        Contig N20: 2516
        Contig N30: 1841
        Contig N40: 1315
        Contig N50: 895 #at least half of the transcripts are this size in bp

        Median contig length: 357
        Average contig: 624.22
        Total assembled bases: 332798719

#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3422
        Contig N20: 2229
        Contig N30: 1487
        Contig N40: 931
        Contig N50: 625

        Median contig length: 316
        Average contig: 525.95
        Total assembled bases: 154171724

#The N10 through N50 values are shown computed based on all assembled contigs. In this example, 10% of the assembled bases are found in transcript contigs at least 3,669 bases in length (N10 value), and the N50 value indicates that at least half the assembled bases are found in contigs that are at least 895 bases in length.

#We use metazoa

##################### 6.8) BUSCO ###########################

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=crubrum_busco
#SBATCH --mem=100G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=30            # default = 1, max = 48
#SBATCH --output=busco-crubrum.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=busco-crubrum.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen

####### BUSCO module already loaded and included in my PATH #######

module load hmmer/3.3

########### Do interesting things here:

busco -i /home/sramirezc/sramirezc/crubrum_rnaseq/c.rubrum_reference/crubrum.trinity.fasta -o /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/BUSCO -l ./busco_downloads/lineages/metazoa_odb10 -f -m tran -c 30 --metaeuk_parameters="--disk-space-limit=500M,--remove-tmp-files=1" --metaeuk_rerun_parameters="--disk-space-limit=500M,--remove-tmp-files=1"

####### Script ends here #####

# -f --> force the analysis despite there are outpufiles already
# --update-data --> update the linage version automatically

## output

2022-12-13 14:46:55 INFO:

        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:99.3%[S:17.3%,D:82.0%],F:0.2%,M:0.5%,n:954     |
        |947    Complete BUSCOs (C)                       |
        |165    Complete and single-copy BUSCOs (S)       |
        |782    Complete and duplicated BUSCOs (D)        |
        |2      Fragmented BUSCOs (F)                     |
        |5      Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
        
2022-12-13 14:46:55 INFO:       BUSCO analysis done. Total running time: 1370 seconds
2022-12-13 14:46:55 INFO:       Results written in /mnt/lustre/bio/users/sramirezc/crubrum_rnaseq/assembl_quality/BUSCO/home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/BUSCO

#provide a quantitative assessment of the completeness in terms of expected gene content of a transcriptome, 
#a high level of duplication may be explained by a recent whole duplication event (biological) or a chimeric assembly of haplotypes (technical). Transcriptomes and protein sets that are not filtered for isoforms will lead to a high proportion of duplicates. Therefore you should filter them before a BUSCO analysis. 

############################ Further annotation to merge with differential expression analysis #########################

#our final reference
crubrum_final.transdecoder.fasta

#N. of transcripts: 
grep -c ">" crubrum_final.transdecoder.fasta
#124,373

#Script used for Swissprot data base

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=blastx_swiss
#SBATCH --mem=500G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=48            # default = 1, max = 48
#SBATCH --output=blastx_swissprot.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_swissprot.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen


###### load modules

module load blast/2.7.1

# do the interesting things here:

blastx -query crubrum_final.transdecoder.fasta -db swissprot.fa -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/blast/swissprot_res/crubrum_blastx_swissprot.outfmt14.txt -max_hsps 1 -evalue 1e-5 -num_threads 48 -max_target_seqs 5 -outfmt 14

#We use blastx since we have a query with translated protenis obtained after filtering with transdecoder. This way is more accurate and exact than with nucleotides

### Now with softcorals from ncbi

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=blastx_softcoral
#SBATCH --mem=500G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=48            # default = 1, max = 48
#SBATCH --output=blastx_softcoral.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_softcoral.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen


###### load modules

module load blast/2.7.1

# do the interesting things here:

blastx -query crubrum_final.transdecoder_def.fasta -db softcoralsdb -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/databases/softcorals/result/c.rubrum_blastx_uniprot.outfmt14.txt -max_hsps 1 -evalue 1e-5 -num_threads 48 -max_target_seqs 5 -outfmt 14

#finally with corals from ensembl:

#!/bin/bash
#SBATCH --account=mbc1           #bank_account
#SBATCH --job-name=blastx_coral
#SBATCH --mem=50G                 # Memory in MB put the G
#SBATCH --ntasks=1                    # related with n of nodes
#SBATCH --cpus-per-task=48            # default = 1, max = 48
#SBATCH --output=blastx_coral.log.out     # File to which standard out will be written | path for output files
#SBATCH --error=blastx_coral.err      # File to which standard err will be written | provide path
##SBATCH --array=1-72%4                  #se pone un rango con el n de elementos que se tienen


###### load modules

module load blast/2.7.1

# do the interesting things here:

blastx -query crubrum_final.transdecoder_def.fasta -db coralsdb -out /home/sramirezc/sramirezc/crubrum_rnaseq/assembl_quality/coral_embl/result_coralblastx_fmt14/c.rubrum_blastx_coral.outfmt14.txt -max_hsps 1 -evalue 1e-5 -num_threads 48 -max_target_seqs 5 -outfmt 14

#The results will be uploaded to Omicsbox to assign GO terms and interpro IDs.
