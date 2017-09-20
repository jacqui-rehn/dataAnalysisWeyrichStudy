#!/bin/bash

#count *bwa_bam and *rmdup.bam and *split.bam files at each MAPQ after length filtering

#Specify variables
ROOTDIR=/home/a1698312
MAPDIR=$ROOTDIR/weyrich/filteredMapData
LQ_MAPDIR=$MAPDIR/lowQualMapData
HQ_MAPDIR=$MAPDIR/highQualMapData

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

#Count reads in *bwa.bam at different MAPQ scores
for bwa_file in *bwa.bam
  do
    echo "Counting reads in ${bwa_file}"
    MAPCOUNT=$(samtools view ${bwa_file} | cut -f5 | sort | uniq -c)
    echo -e "#Read\tMAPQ" > ${bwa_file%%.bam}_count.txt
    echo -e "${MAPCOUNT}" >> ${bwa_file%%.bam}_count.txt
  done
  
#Count reads in rmdup.bam at different MAPQ scores
for rmdup_file in *rmdup.bam
do
  echo "Counting reads in ${rmdup_file}"
  MAPCOUNT2=$(samtools view ${rmdup_file} | cut -f5 | sort | uniq -c)
  echo -e "#Read\tMAPQ" > ${rmdup_file%%.bam}_count.txt
  echo -e "${MAPCOUNT2}" >> ${rmdup_file%%.bam}_count.txt
done


#Generate text file for storing alignment count data for each genome
if [ ! -f filtered_split_count.txt ]
then
  echo -e 'Creating file filtered_split_count.txt'
  echo -e "#Reads  MAPQ" > filtered_split_count.txt
else
  echo  'filtered_split_count.txt already exists'
fi

#Count reads at different MAPQ scores for each rmdup.bam
for split_file in *split.bam
  do
    echo "Counting reads in ${split_file}"
    MAPCOUNT3=$(samtools view ${split_file} | cut -f5 | sort | uniq -c)
    echo -e "${split_file}" >> filtered_split_count.txt
    echo -e "${MAPCOUNT3}" >> filtered_split_count.txt
  done
  
  



