#Alignment Count Data

#load packages
library(dplyr)
library(readr)
library(magrittr)
library(tibble)
library(stringr)
library(reshape2)
library(ggplot2)
library(data.table)

#Read in txt file with fastq count data and assign to object
fastqCount <- read_delim("trimData/collapsed_read_count.txt", delim = "\t", 
           skip = 1, col_names = FALSE) %>%
  set_colnames(c("sampleID", "mergedReads"))

#Edit file name
fastqCount$sampleID <- gsub('2NoAdapt_', '', fastqCount$sampleID)
fastqCount$sampleID <- gsub('_R1R2_Collapsed', '', fastqCount$sampleID)
fastqCount$sampleID <- gsub('L7_lTACTG_rCTCGA', '', fastqCount$sampleID)
fastqCount$sampleID <- gsub('L7_lAAGAG_rNONE', '', fastqCount$sampleID)

#Read in text file with alingment count data and assign to object
alnCount <- read_delim("alnData/aligned_read_count.txt", delim = "\t", skip = 1, 
                       col_names = FALSE) %>% 
  set_colnames(c("fileName", "alnReads", "alnQ10", "alnQ20", "alnQ30"))

#Edit fileNames
alnCount$fileName <- gsub('2NoAdapt_', '', alnCount$fileName)
alnCount$fileName <- gsub('L7_lTACTG_rCTCGA', '', alnCount$fileName)
alnCount$fileName <- gsub('L7_lAAGAG_rNONE', '', alnCount$fileName)
alnCount$fileName <- gsub('_split', '', alnCount$fileName)

#Reshape data to separate sampleID and genome aligning with from fileName information
alnCount <- alnCount %>% mutate(genome = str_extract(fileName, 
            "(aoris|bwa|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|rmdup|smutans|tdenticola|tforsythia|moralis)"), 
            sampleID = str_replace(fileName, 
            "_(aoris|bwa|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|rmdup|smutans|tdenticola|tforsythia|moralis)", ""))
alnCount <- alnCount[, c(7,6,2:5)]
  
#extract information about the total number of reads aligned
bwaCount <- alnCount %>% filter(genome == "bwa") %>% select(-genome, -alnQ10, -alnQ20, -alnQ30)
#extract information about the total number of reads aligning after rmdup
rmdupCount <- alnCount %>% filter(genome == "rmdup") %>% select(-genome)
#Combine fastqCount, bwaCount and rmdupCount into one data frame
totalCount <- fastqCount %>% left_join(bwaCount) %>% left_join(rmdupCount, by = "sampleID")
#Edit column headings
totalCountHeadings <- c("sampleID", "totalMerged", "totalAln", "totalAlnRmdup", 
                        "alnRmdupQ10", "alnRmdupQ20", "alnRmdupQ30")
names(totalCount) <- totalCountHeadings
#Calculate proportion of reads aligned + proportion of duplicates removed and add to table
totalCount <- totalCount %>% mutate(propAln = (totalAln/totalMerged)*100) %>% 
  mutate(propDupRm = ((totalAln-totalAlnRmdup)/totalAln)*100)
totalCount <- totalCount[,c(1:3,8,9,4:7)]
#Calculate proprotion of aligned reads at different quality scores and add to table
totalCount <- totalCount %>% mutate(propQ10 = (alnRmdupQ10/totalAlnRmdup)*100) %>% 
  mutate(propQ20 = (alnRmdupQ20/totalAlnRmdup)*100) %>% 
  mutate(propQ30 = (alnRmdupQ30/totalAlnRmdup)*100)
totalCount <- totalCount[, c(1:7,10,8,11,9,12)]

#Remove bwa data from alnCount data frame
alnCount <- alnCount %>% filter(genome != "bwa")
#Edit sampleID names so that all are the same
alnCount$sampleID <- gsub('Elsidron1', 'ELSIDRON1', alnCount$sampleID)
alnCount$sampleID <- gsub('ModernL7', 'Modern', alnCount$sampleID)
alnCount <- alnCount[with(alnCount, order(sampleID)), ]

#Calculate proportion of reads aligning to each genome using totalAlnRmdup
##ELSIDRON1 totalAlnRmdup = 547,017
##Modern totalAlnRmdup = 582,311
alnCount <- rmdupCount %>% select(-alnQ10, -alnQ20, -alnQ30) %>% left_join(alnCount, by = "sampleID")
#edit heading names
alnCountHeadings <- c("sampleID", "alnReads", "genome", "alnQ0", "alnQ10", "alnQ20", "alnQ30")
names(alnCount) <- alnCountHeadings
#Convert rmdup to all
alnCount$genome <- gsub('rmdup', 'all', alnCount$genome)
#sort columns so that genomes listed in alphabetical order
alnCount <- alnCount[with(alnCount, order(sampleID, genome)), ]

#Convert alnCount from wide to long format 
alnCountLong <- alnCount %>% select(-alnReads) %>% 
  melt(id.vars = c("sampleID", "genome"), variable.name = "qual.cutoff", value.name = "count")
#substitute alnQ0 for a value
alnCountLong$qual.cutoff <- gsub('alnQ', '', alnCountLong$qual.cutoff)
#Convert qual.cutoff from chr to as.numeric
alnCountLong$qual.cutoff <- as.numeric(alnCountLong$qual.cutoff)

#plot line graph representing number aligning 
ElsidronCountPlot1 <- alnCountLong %>% filter(sampleID == "ELSIDRON1") %>% 
  ggplot(aes(x=qual.cutoff, y=count, group=genome, colour=genome)) + geom_line() +
  theme_bw() + xlab("Aligment quality score cutoff") +
  scale_y_continuous(labels = scales::comma) +
  scale_colour_manual(values=c("firebrick1", "slateblue3", "darkseagreen", "goldenrod3", "plum3",
                               "coral2", "deepskyblue4", "brown1", "seagreen3", "violetred3", 
                               "royalblue3")) +
  ylab("Number of reads aligning") + 
  ggtitle("Effect of quality score on the number of\naligned reads in Elsidron1 sample")

#same plot as above but with limits set to ensure colours stay consistent
alnCountLong %>% filter(sampleID == "ELSIDRON1") %>% 
  ggplot(aes(x=qual.cutoff, y=count, group=genome, colour=genome)) + geom_line() +
  theme_bw() + xlab("Aligment quality score cutoff") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_discrete(drop = TRUE, limits = levels(alnCountLong$genome)) +
  ggtitle("Effect of quality score on the number of\naligned reads in Elsidron1 sample")

#Test to see if colours consistent when plot with Modern sample
alnCountLong %>% filter(sampleID == "Modern") %>% 
  ggplot(aes(x=qual.cutoff, y=count, group=genome, colour=genome)) + geom_line() +
  theme_bw() + xlab("Aligment quality score cutoff") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_discrete(drop = TRUE, limits = levels(alnCountLong$genome)) +
  ggtitle("Effect of quality score on the number of\naligned reads in Modern sample")

#Test to see if I can use facet-wrap to create simultaneous plots
alnCountLong %>% ggplot(aes(x=qual.cutoff, y=count, group=genome, colour=genome)) + geom_line() +
  theme_bw() + xlab("Aligment quality score cutoff") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_discrete(drop = TRUE, limits = levels(alnCountLong$genome)) +
  facet_wrap(~sampleID)

#Subsetted plot focusing on genomes with large proportion of reads aliging to
alnCountLong %>% filter(sampleID == "ELSIDRON1") %>% filter(count > 75000) %>%
  ggplot(aes(x=qual.cutoff, y=count, group=genome, colour=genome)) + geom_line() +
  theme_bw() + xlab("Aligment quality score cutoff") +
  scale_y_continuous(labels = scales::comma) +
  ylab("Number of reads aligning") + 
  ggtitle("Effect of quality score on the number of\naligned reads in Elsidron1 sample")

#Subset with genomes to which a moderate number of reads align
alnCountLong %>% filter(sampleID == "ELSIDRON1") %>% subset(count>8000 & count<75000) %>%
  ggplot(aes(x=qual.cutoff, y=count, group=genome, colour = genome)) + geom_line() +
  theme_bw() + xlab("Aligment quality score cutoff") +
  scale_y_continuous(labels = scales::comma) +
  ylab("Number of reads aligning") + 
  ggtitle("Effect of quality score on the number of\naligned reads in Elsidron1 sample") 

#Subset of genomes to which a small number of reads align
alnCountLong %>% filter(sampleID == "ELSIDRON1") %>% filter(count < 10000) %>%
  ggplot(aes(x=qual.cutoff, y=count, group=genome, colour=genome)) + geom_point() +
  stat_smooth(method=lm, se=FALSE) + theme_bw() + xlab("Aligment quality score cutoff") +
  ylab("Number of reads aligning") + 
  ggtitle("Effect of quality score on the number of\naligned reads in Elsidron1 sample")

#Convert alnCounts to alnProps and assign to an object
alnProps <- alnCount %>% mutate(alnQ0 = (alnQ0/alnReads)*100) %>% 
  mutate(alnQ10 = (alnQ10/alnReads)*100) %>% 
  mutate(alnQ20 = (alnQ20/alnReads)*100) %>% 
  mutate(alnQ30 = (alnQ30/alnReads)*100)
#Convert from wide to long format
alnPropsLong <- alnProps %>% select(-alnReads) %>% 
  melt(id.vars = c("sampleID", "genome"), variable.name = "qual.cutoff", 
       value.name = "prop", variable.facto = FALSE)

#Convert qual.cutoff from factor to character
alnPropsLong$qual.cutoff <- as.character(alnPropsLong$qual.cutoff)
#Remove alnQ from the character vector qual.cutoff
alnPropsLong$qual.cutoff <- gsub('alnQ', '', alnPropsLong$qual.cutoff)
#Convert qual.cutoff from character to numeric
alnPropsLong$qual.cutoff <- as.numeric(alnPropsLong$qual.cutoff)

#Plot proportions alinging for each cutoff for each genome
alnPropsLong %>% filter(sampleID == "ELSIDRON1") %>% 
  ggplot(aes(x=qual.cutoff, y=prop, group=genome, colour=genome)) + 
  geom_line() + theme_bw() + xlab("Aligment quality score cutoff") +
  ylab("% of total aligned reads") + 
  ggtitle("Effect of quality score on proportion of Elsidron1\nreads aligning to each genome")

alnPropsLong %>% filter(sampleID == "Modern") %>%
  ggplot(aes(x=qual.cutoff, y=prop, group=genome, colour=genome)) + 
  geom_line() + theme_bw() + xlab("Aligment quality score cutoff") +
  ylab("% of total aligned reads") + 
  ggtitle("Effect of quality score on proportion of Modern\nreads aligning to each genome")

#Plot stacked bar plot indicating number of reads aligning to each genome
alnCountLong %>% filter(genome != "all") %>% ggplot(aes(x=qual.cutoff, y=count, fill=genome)) + 
  geom_bar(stat = "identity") + facet_wrap(~sampleID)

