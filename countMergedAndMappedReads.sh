#!/bin/bash

#count merged_fastq and _rmdup.bam and _split.bam files

#Specify variables
ROOTDIR=/home/a1698312
TRIMDIR=$ROOTDIR/weyrich/trimData
#ALNDIR=$ROOTDIR/weyrich/alnData
#QUALALNDIR=$ROOTDIR/weyrich/highQualAlnData
MAPDIR=$ROOTDIR/weyrich/alnData
LOW_Q_DIR=$MAPDIR/lowQualMapData
HIGH_Q_DIR=$MAPDIR/highQualMapData

########### Count merged Reads ###############

#Change into directory where Collapsed_fastq files located
#if [ -d ${TRIMDIR} ]
#then
#  echo "Changing to trimData directory"
#  cd ${TRIMDIR}
#else
#  echo "Cannot find ${TRIMDIR}"
#exit1
#fi

#Generate text file for storing merged count data
#if [ ! -f fastq_read_count.txt ]
#then
#  echo -e 'Creating file fastq_read_count.txt'
#  echo -e 'fileName\t#Reads' > fastq_read_count.txt
#else
#  echo  'fastq count file already exists'
#fi

#Count merged reads in fastq files and print to text file
#for Collapsed_fastq in *_Collapsed.fastq.gz
#  do
#    echo "Counting number of merged reads in ${Collapsed_fastq}"
#    MERGECOUNT=$(zcat ${Collapsed_fastq} | egrep -c '^@M_HWI')
#    echo -e "${Collapsed_fastq%%_R1R2_Collapsed.fastq.gz}\t${MERGECOUNT}" >> fastq_read_count.txt
#  done

########## Count aligned reads ############

#Change into directory where aln_files located
if [ -d ${MAPDIR} ]
then
  echo "Changing to mapData directory"
  cd ${MAPDIR}
else
  echo "Cannot find ${MAPDIR}"
exit1
fi

#Generate text file for storing alignment count data
if [ ! -f bam_read_count.txt ]
then
  echo -e 'Creating file bam_read_count.txt'
  echo -e "readCount\tMAPQ" > bam_read_count.txt
else
  echo  'Bam count file already exists'
fi

#Count total number of reads aligned at different quality scores
#for bam_file in *.bam
#  do
#    echo "Counting reads in ${bam_file}"
#    ALNCOUNT=$(samtools view -c ${bam_file})
#    ALNCOUNTQ10=$(samtools view -q 10 -c ${bam_file})
#    ALNCOUNTQ20=$(samtools view -q 20 -c ${bam_file})
#    ALNCOUNTQ30=$(samtools view -q 30 -c ${bam_file})
#    MAPCOUNT=$(samtools view ${bam_file} | cut -f5 | sort | uniq -c)
#    echo -e "${bam_file%%.bam}\t${ALNCOUNT}\t${ALNCOUNTQ10}\t${ALNCOUNTQ20}\t${ALNCOUNTQ30}" >> aligned_read_count.txt
#    echo -e "${bam_file}" >> bam_read_count.txt
#    echo -e "${MAPCOUNT}" >> bam_read_count.txt
#  done

#Count reads at different MAPQ scores for each rmdup.bam
for bam_file in *bwa.bam
  do
    echo "Counting reads in ${bam_file}"
    MAPCOUNT=$(samtools view ${bam_file} | cut -f5 | sort | uniq -c)
    echo -e "${bam_file}" >> bam_read_count.txt
    echo -e "${MAPCOUNT}" >> bam_read_count.txt
  done

#Generate text file for storing alignment count data
if [ ! -f rmdup_read_count.txt ]
then
  echo -e 'Creating file rmdup_read_count.txt'
  echo -e "readCount\tMAPQ" > rmdup_read_count.txt
else
  echo  'rmdup count file already exists'
fi

#Count reads at different MAPQ scores for each rmdup.bam
for rmdup_file in *rmdup.bam
  do
    echo "Counting reads in ${rmdup_file}"
    MAPCOUNT=$(samtools view ${rmdup_file} | cut -f5 | sort | uniq -c)
    echo -e "${rmdup_file}" >> rmdup_read_count.txt
    echo -e "${MAPCOUNT}" >> rmdup_read_count.txt
  done

##############Count aligned reads in split.bam files################

#Generate text file for storing split alignment count data
#if [ ! -f split_aligned_read_count.txt ]
#then
#  echo -e 'Creating file split_aligned_read_count.txt'
#  echo -e "fileName\t#alnReads" > split_aligned_read_count.txt
#else
#  echo  'split_aligned_read_count.txt already exists'
#fi

#Count total number of reads aligned in each split file with no quality filtering
#for split_file in *split.bam
#  do
#    echo "Counting reads in ${split_file}"
#    SPLITCOUNT=$(samtools view -c ${split_file})
#    echo -e "${split_file%%_split.bam}\t${SPLITCOUNT}" >> split_aligned_read_count.txt
#  done

#Change into directory where highQualAln_files located
#if [ -d ${QUALALNDIR} ]
#then
#  echo "Changing to highQualAlnData directory"
#  cd ${QUALALNDIR}
#else
#  echo "Cannot find ${QUALALNDIR}"
#exit1
#fi

#Generate text file for storing split alignment count data
#if [ ! -f highQual_split_aligned_read_count.txt ]
#then
#  echo -e 'Creating file highQual_split_aligned_read_count.txt'
#  echo -e "fileName\t#alnReads" > highQual_split_aligned_read_count.txt
#else
#  echo  'highQual_split_aligned_read_count.txt already exists'
#fi

#Count total number of reads aligned in each split file with -q 30 filter
#for Q30split_file in *Q30_split.bam
#  do
#    echo "Counting reads in ${Q30split_file}"
#    Q30SPLITCOUNT=$(samtools view -c ${Q30split_file})
#    echo -e "${Q30split_file%%_Q30_split.bam}\t${Q30SPLITCOUNT}" >> highQual_split_aligned_read_count.txt
#  done

#Generate text file for storing split alignment count data
if [ ! -f split_read_count.txt ]
then
  echo -e 'Creating file split_read_count.txt'
  echo -e "readCount\tMAPQ" > split_read_count.txt
else
  echo  'split_read_count.txt already exists'
fi

#Count reads at different MAPQ scores for each split.bam
for split_file in *split.bam
  do
    echo "Counting reads in ${split_file}"
    MAPCOUNT=$(samtools view ${split_file} | cut -f5 | sort | uniq -c)
    echo -e "${split_file}" >> split_read_count.txt
    echo -e "${MAPCOUNT}" >> split_read_count.txt
  done

########### Count split.bam at each MAPQ for lowQualMapData ###################

#Change into directory where low_qual_split_files located
#if [ -d ${LOW_Q_DIR} ]
#then
#  echo "Changing to lowQualMapData directory"
#  cd ${LOW_Q_DIR}
#else
#  echo "Cannot find ${LOW_Q_DIR}"
#exit1
#fi

#Generate text file for storing alignment count data
#if [ ! -f low_qual_split_count.txt ]
#then
#  echo -e 'Creating file low_qual_split_count.txt'
#  echo -e "#Reads  MAPQ" > low_qual_split_count.txt
#else
#  echo  'low_qual_split_count.txt already exists'
#fi

#Count reads at different MAPQ scores for each rmdup.bam
#for low_qual_split_file in *split.bam
#  do
#    echo "Counting reads in ${low_qual_split_file}"
#    MAPCOUNT=$(samtools view ${low_qual_split_file} | cut -f5 | sort | uniq -c)
#    echo -e "${low_qual_split_file}" >> low_qual_split_count.txt
#    echo -e "${MAPCOUNT}" >> low_qual_split_count.txt
#  done

########### Count reads in highQualMapData split.bam ###################

#Change into directory where highQualAln_files located
#if [ -d ${HIGH_Q_DIR} ]
#then
#  echo "Changing to highQualMapData directory"
#  cd ${HIGH_Q_DIR}
#else
#  echo "Cannot find ${HIGH_Q_DIR}"
#exit1
#fi

#Generate text file for storing split alignment count data
#if [ ! -f high_qual_split_count.txt ]
#then
#  echo -e 'Creating file high_qual_split_count.txt'
#  echo -e "fileName\tNo.Reads" > high_qual_split_count.txt
#else
#  echo  'high_qual_split_count.txt already exists'
#fi

#Count total number of reads aligned in each split file with -q 30 filter
#for high_qual_split_file in *split.bam
#  do
#    echo "Counting reads in ${high_qual_split_file}"
#    MAPCOUNT2=$(samtools view -c ${high_qual_split_file})
#    echo -e "${high_qual_split_file%%_split.bam}\t${MAPCOUNT2}" >> high_qual_split_count.txt
#  done