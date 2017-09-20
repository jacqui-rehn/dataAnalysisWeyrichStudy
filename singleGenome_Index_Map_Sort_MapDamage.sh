#!/bin/bash

#USAGE: Requires trimmed_fastq files in specified directory
#       Specify variable for location of Ref Seq genomes for downloading and alignment

#Specify variables
ROOTDIR=/home/a1698312
TRIMDIR=$ROOTDIR/weyrich/trimData
MAPDIR=$ROOTDIR/weyrich/smutans
LOGFILE=$ROOTDIR/weyrich/mapDamageLog.txt
REF=S.mutans
SAMPLE=Elsidron1

##### bwa_build alignment index and align to smutans####

#Add time stamp to logfile
echo -e "$(date -u)\t start - build bwa index and align to smutans" >> ${LOGFILE}

#Change into mapData directory
if [ -d ${MAPDIR} ]
then
  echo "Changing into ${MAPDIR}"
  cd ${MAPDIR}
else
  echo "Cannot find ${MAPDIR}"
exit1
fi

#Download smutans refSeq genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/465/GCF_000007465.2_ASM746v2/GCF_000007465.2_ASM746v2_genomic.fna.gz -O smutans.fna.gz

#unzip fasta files
gunzip *fna.gz

#build-index for alignment
bwa index -p bwaidx smutans.fna

#Change into directory where trimmed fastq files located
if [ -d ${TRIMDIR} ]
then
  echo "Changing to trimData directory"
  cd ${TRIMDIR}
else
  echo "Cannot find ${TRIMDIR}"
exit1
fi

#bwa alignment of collapsed reads
echo "Aligning 2NoAdapt_ELSIDRON1L7_lTACTG_rCTCGA_R1R2_Collapsed.fastq.gz"
bwa aln -n 0.01 -o 2 -l 1024 -t 4 $MAPDIR/bwaidx 2NoAdapt_ELSIDRON1L7_lTACTG_rCTCGA_R1R2_Collapsed.fastq.gz > Elsidron1_S.mutans_MAPPED.sai


#Convert .sai alignment file to bam format with the sam header. Exclude unmapped reads.
echo "Converting Elsidron1_S.mutans_MAPPED.sai to bam format"
    bwa samse $MAPDIR/bwaidx \
              Elsidron1_S.mutans_MAPPED.sai \
              2NoAdapt_ELSIDRON1L7_lTACTG_rCTCGA_R1R2_Collapsed.fastq.gz | \
                samtools view -bSh -F0x4 -> $MAPDIR/Elsidron1_S.mutans_bwa.bam

#Add time stamp to logfile
echo -e "$(date -u)\t finish - build bwa index and align to smutans" >> ${LOGFILE}

################### sambamba sort and rmdup and count #####################

#Add time stamp to logfile
echo -e "$(date -u)\t start - sambamba sort/rmdup and count" >> ${LOGFILE}

#Change into mapData directory
if [ -d ${MAPDIR} ]
then
  echo "Changing to mapData directory"
  cd ${MAPDIR}
else
  echo "Cannot find ${MAPDIR}"
exit1
fi

#sort bam file
echo "Sorting bam file for ${bam_file}"
sambamba sort -o Elsidron1_S.mutans_sorted.bam Elsidron1_S.mutans_bwa.bam

#remove duplicates
echo "Removing duplicates ${sort_file}"
sambamba markdup -r Elsidron1_S.mutans_sorted.bam Elsidron1_S.mutans_rmdup.bam

#Remove _sorted.bam.bai files as no longer needed
rm *_sorted.bam.bai

#Generate text file for count data
if [ ! -f read_count.txt ]
then
  echo -e 'Creating file read_count.txt'
  echo -e "readCount\tMAPQ" > read_count.txt
else
  echo  'read count file already exists'
fi

#Count reads at different MAPQ scores for each *.bam
for bam_file in *.bam
  do
    echo "Counting reads in ${bam_file}"
    MAPCOUNT=$(samtools view ${bam_file} | cut -f5 | sort | uniq -c)
    echo -e "${bam_file}" >> read_count.txt
    echo -e "${MAPCOUNT}" >> read_count.txt
  done

echo -e "$(date -u)\t finish - sambamba sort/rmdup and count" >> ${LOGFILE}

################# mapDamage #########################

echo -e "$(date -u)\t start - mapDamage Elsidron1_S.mutans" >> ${LOGFILE}

echo "Running mapDamage on Elsidron1_S.mutans_rmdup.bam"
mapDamage -i Elsidron1_S.mutans_rmdup.bam -r $MAPDIR/smutans.fna

echo -e "$(date -u)\t finish - mapDamage Elsidron1_S.mutans" >> ${LOGFILE}
