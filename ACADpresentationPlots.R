#load packages
library(dplyr)
library(readr)
library(magrittr)
library(tibble)
library(stringr)
library(reshape2)
library(ggplot2)
library(data.table)
library(scales)
library(lme4)
library(lmerTest)

############## Original Data - 2 Samples; 10 genomes ####################

### Count Data table

#Read in txt file with fastq count data and assign to object
fastqCount <- read_delim("trimData/collapsed_read_count.txt", delim = "\t", 
                         skip = 1, col_names = FALSE) %>%
  set_colnames(c("fileName", "fastqCount"))

#Edit file name to include only sampleID
fastqCount <- colsplit(fastqCount$fileName, "_", names=c("adapters", "sampleID", "extra")) %>% 
  bind_cols(fastqCount) %>% 
  select(-fileName, -adapters, -extra)

#Read in bamCount csv file
bwaCount <- read.csv(file="alnData/bam_read_count.txt", sep="", skip = 1, header = FALSE, col.names = c("bwaCount", "MAPQ"))
#Split at .bam to generate a list
bwaCount <- bwaCount %>% mutate(bam = grepl("bam", bwaCount), fileNo = cumsum(bam)) %>% split(f = .$fileNo)
#Counts listed in same order files are listed within directory
#Therefore, generate a list of all bamFiles that were counted
bwaFiles <- list.files("alnData/", pattern = "_bwa.bam$", full.names = FALSE)
#assign this list as names of files in bamCount
names(bwaCount) <- bwaFiles
#Bind_rows of list, taking list names and re-inserting as fileName, then remove unnecessary columns and rows
bwaCount <- bwaCount %>% bind_rows(.id = "fileName") %>% select(-bam, -fileNo) %>% filter(MAPQ != "NA")
#split fileName into adapters, sampleID and additional info
bwaCount <- colsplit(bwaCount$fileName, "_", names=c("adapters", "sampleID", "extra")) %>% 
  bind_cols(bwaCount) %>% 
  select(-fileName, -adapters, -extra)
#Convert bwaCount variable from factor to numeric
bwaCount %>% mutate_if(is.factor, as.character) -> bwaCount
bwaCount$bwaCount <- as.numeric(bwaCount$bwaCount)

#Read in rmdupCount csv file
rmdupCount <- read.csv(file="alnData/rmdup_read_count.txt", sep="", skip = 1, header = FALSE, col.names = c("rmdupCount", "MAPQ"))
#Split at .bam to generate a list
rmdupCount <- rmdupCount %>% mutate(bam = grepl("bam", rmdupCount), fileNo = cumsum(bam)) %>% split(f = .$fileNo)
#Counts listed in same order files are listed within directory
#Therefore, generate a list of all bamFiles that were counted
rmdupFiles <- list.files("alnData/", pattern = "_rmdup.bam$", full.names = FALSE)
#assign this list as names of files in bamCount
names(rmdupCount) <- rmdupFiles
#Bind_rows of list, taking list names and re-inserting as fileName, then remove unnecessary information
rmdupCount <- rmdupCount %>% bind_rows(.id = "fileName") %>% select(-bam, -fileNo) %>% filter(MAPQ != "NA")
#split fileName into adapters, sampleID and additional info
rmdupCount <- colsplit(rmdupCount$fileName, "_", names=c("adapters", "sampleID", "extra")) %>% 
  bind_cols(rmdupCount) %>% 
  select(-fileName, -adapters, -extra)
#Convert bwaCount variable from factor to numeric
rmdupCount %>% mutate_if(is.factor, as.character) -> rmdupCount
rmdupCount$rmdupCount <- as.numeric(rmdupCount$rmdupCount)

######Create table summarising total fastq, bwa.bam and rmdup.bam for each sample#####
rmdupCount %>% select(-MAPQ) %>% group_by(sampleID) %>% summarise_each(funs(sum)) -> totalRmdupCount
bwaCount %>% select(-MAPQ) %>% group_by(sampleID) %>% summarise_each(funs(sum)) -> totalBwaCount
totalCount <- left_join(fastqCount, totalBwaCount, by = "sampleID")
totalCount <- left_join(totalCount, totalRmdupCount, by = "sampleID")
##From this calculate the total number of duplicate reads identified in each sample and add to totalCount
totalCount %>% mutate(dupCount = bwaCount - rmdupCount) -> totalCount
totalCount <- totalCount[, c(1:3,5,4)]

#create list of sample names for each ID
sample_names <- c(
  `ELSIDRON1L7` = "Elsidron 1",
  `ModernL7` = "Modern"
)

###Plot number of reads sequenced and aligned for each sample
countPlot1 <- totalCount %>% 
  select(sampleID, fastqCount, bwaCount) %>% 
  melt(id.vars = c("sampleID"), variable.name = "counting", value.name = "count") %>% 
  ggplot(aes(x="", y=count, fill=counting)) + 
  geom_bar(width = 1, stat = "identity") + 
  scale_y_continuous(labels = scales::comma) + 
  theme_bw() +
  guides(fill=guide_legend(title=NULL)) + 
  scale_fill_manual(values = c("#3399FF", "#FF6666"), 
                    breaks=c("bwaCount", "fastqCount"), 
                    labels=c("Aligned", "Not aligned")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank()) + 
  facet_wrap(~sampleID, labeller = as_labeller(sample_names))

#Create a blank theme to be applied to pie charts
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="plain", hjust = 0.5)
  )

#Calculate % of reads aligned/unaligned, duplicate/nonduplicate and add to totalCount
totalCount <- totalCount %>% 
  mutate(propAln = (bwaCount/fastqCount)*100, 
         propUnAln = ((fastqCount-bwaCount)/fastqCount)*100, 
         propDup = (dupCount/bwaCount)*100, 
         propNonDup = (rmdupCount/bwaCount)*100)
#round % to 2 decimal places
totalCount[,6:9] <- round(totalCount[,6:9],2)

countPlot2 <- totalCount %>% 
  select(sampleID, propAln, propUnAln) %>% 
  melt(id.vars = "sampleID") %>% 
  ggplot(aes(x="", y=value, fill=variable)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) + 
  blank_theme +
  scale_fill_manual(values = c("#FF6666", "#3399FF"), 
                    breaks=c("bwaCount", "fastqCount"), 
                    labels=c("Aligned", "Not aligned")) +
  theme(axis.text.x = element_blank(), strip.text = element_text(size = 11)) + 
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  facet_wrap(~sampleID, ncol = 1, labeller = as_labeller(sample_names))

#Put both plots side by side
#pdf("plotsACADpres/initialData_AlnCounts.pdf")
grid.arrange(countPlot1, countPlot2, ncol=2, widths=c(1,0.5))
#dev.off()

#Plot the proportion of duplicates removed from each sample
#pdf("plotsACADpres/initialData_DupCounts.pdf")
totalCount %>%
  select(sampleID, propDup, propNonDup) %>% 
  melt(id.vars = "sampleID") %>% 
  ggplot(aes(x="", y=value, fill=variable)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) + 
  blank_theme +
  scale_fill_manual(values = c("#FF6666", "#3399FF"),
                    name = "",
                    breaks=c("propDup", "propNonDup"), 
                    labels=c("Duplicate", "Non-duplicate")) +
  theme(axis.text.x = element_blank()) + 
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  facet_wrap(~sampleID, labeller = as_labeller(sample_names)) + 
  theme(strip.text = element_text(size = 12))
#dev.off()

###plot No Reads by MAPQ, before and after de-duplication
#combine bwaCount and rmdupCount data
MAPQcount <- left_join(bwaCount, rmdupCount)
#calculate #duplicates present at each MAPQ
MAPQcount <- MAPQcount %>% mutate(dupCount = bwaCount - rmdupCount)
#re-order cols
MAPQcount <- MAPQcount[, c(1,3,2,4,5)]
#Collate counts into bins
MAPQcount$MAPQrange <- cut(MAPQcount$MAPQ, breaks = c(0,10,20,30,40), 
                           labels = c("0-9", "10-19", "20-29", "30-40"), 
                           right = FALSE)
#split data frame by sampleID and summarise counts within each MAPQrange and return to single data frame
MAPQcount <- MAPQcount %>% 
  split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
      select(bwaCount, rmdupCount, dupCount, MAPQrange) %>% 
      group_by(MAPQrange) %>% 
      summarise_each(funs(sum))}) %>% 
  bind_rows(.id = "sampleID")

#Plot MAPQ for each value
pdf("plotsACADpres/initialData_ReadsByMAPQ.pdf")
MAPQcount %>%
melt(id.vars = c("sampleID", "MAPQrange")) %>% 
  ggplot(aes(x=MAPQrange, y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(~sampleID, labeller = as_labeller(sample_names)) + 
  theme_bw() + 
  scale_y_continuous(labels = scales::comma) +
  ylab("") + 
  xlab("MAPQ score") + 
  ggtitle("Number of reads by MAPQ score") + 
  guides(fill=guide_legend(title = NULL)) + 
  scale_fill_manual(values = c("#FF6666", "#66CC66", "#3399FF"), 
                    labels=c("Aligned", "Duplicate", "Remaining"))
dev.off()

###plot prop reads for each MAPQ
#calculate proportions of aln, dup, non-dup for each MAPQrange
MAPQprop <- totalCount %>% 
  select(sampleID, bwaCount, dupCount, rmdupCount) %>% 
  left_join(MAPQcount, by = "sampleID") %>% 
  mutate(propAln = bwaCount.y/bwaCount.x, 
         propDup = dupCount.y/dupCount.x, 
         propRemain = rmdupCount.y/rmdupCount.x) %>% 
  select(1,5,9:11)
#plot proportions
pdf("plotsACADpres/initialData_PropReadsByMAPQ.pdf")
MAPQprop %>% 
  melt(id.vars = c("sampleID", "MAPQrange")) %>% 
  ggplot(aes(x=MAPQrange, y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(~sampleID, labeller = as_labeller(sample_names)) + 
  theme_bw() +
  ylab("") + 
  xlab("MAPQ score") + 
  ggtitle("Proportion of reads by MAPQ score") + 
  guides(fill=guide_legend(title = NULL)) + 
  scale_fill_manual(values = c("#FF6666", "#66CC66", "#3399FF"), 
                    labels=c("Aligned", "Duplicate", "Remaining"))
dev.off()

###plot splitCounts

#Read in splitCount csv file
splitCount <- read.csv(file="alnData/split_read_count.txt", sep="", skip = 1, header = FALSE, col.names = c("splitCount", "MAPQ"))
#Split at .bam to generate a list
splitCount <- splitCount %>% mutate(bam = grepl("bam", splitCount), fileNo = cumsum(bam)) %>% split(f = .$fileNo)
#Counts listed in same order files are listed within directory
#Therefore, generate a list of all bamFiles that were counted
splitFiles <- list.files("alnData/", pattern = "_split.bam$", full.names = FALSE)
#assign this list as names of files in bamCount
names(splitCount) <- splitFiles
#Bind_rows of list, taking list names and re-inserting as fileName, then remove unnecessary columns and rows
splitCount <- splitCount %>% bind_rows(.id = "fileName") %>% select(-bam, -fileNo) %>% filter(MAPQ != "NA")

#edit fileNames to remove adapter/BC information only sampleID and genome
splitCount$fileName <- gsub('2NoAdapt_', '', splitCount$fileName)
splitCount$fileName <- gsub('lTACTG_rCTCGA_', '', splitCount$fileName)
splitCount$fileName <- gsub('lAAGAG_rNONE_', '', splitCount$fileName)

#split fileName into sampleID, genome and bam
splitCount <- colsplit(splitCount$fileName, "_", names=c("sampleID", "genome", "bam")) %>% 
  bind_cols(splitCount) %>% 
  select(-fileName, -bam)
#Convert bwaCount variable from factor to numeric
splitCount %>% mutate_if(is.factor, as.character) -> splitCount
splitCount$splitCount <- as.numeric(splitCount$splitCount)

#edit sampleID to be identical
splitCount$sampleID <- gsub('Elsidron1', 'ELSIDRON1L7', splitCount$sampleID)

#Collate number of reads for each genome in each sample and assign to object
totalSplitCount <- splitCount %>% 
  split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
      select(-sampleID, -MAPQ) %>% 
      group_by(genome) %>% 
      summarise_each(funs(sum))}) %>% 
  bind_rows(.id = "sampleID")

palette <- c("#FF3333", "#3333FF", "#009900", "#FF9900", "#990099", 
             "#33CCCC", "#66CC66", "#FFCC66", "#FF99CC", "#3399FF", 
             "#FF6666", "#9966FF")

#Combine the above object with totalRmdupCount and calculate prop of total each genome represents & plot as pie chart
pdf("plotsACADpres/initialData_relativeAbundance.pdf")
totalSplitCount %>% 
  left_join(totalRmdupCount, by = "sampleID") %>% 
  mutate(prop = splitCount/rmdupCount) %>% 
  ggplot(aes(x="", y=prop, fill=genome)) + 
  geom_bar(colour = "white", stat = "identity") + 
  coord_polar("y", start = 0) + 
  blank_theme + 
  scale_fill_manual(values = palette, 
                    name = "", 
                    labels=c(genome_names)) +
  facet_wrap(~sampleID, labeller = as_labeller(sample_names)) + 
  theme(axis.text.x = element_blank(), strip.text = element_text(size = 12))
dev.off()

#Plot above, but after quality filter (q > 30)
  #Requires that you first determine total rmdup with q >30 in order to calculate the proportions.
splitCountQ30 <- splitCount %>% filter(MAPQ == "37")
rmdupCountQ30 <- rmdupMAPQcount %>% filter(MAPQrange == "30-40") %>% select(-MAPQrange)






#Collate counts into bins
splitCount$MAPQrange <- cut(splitCount$MAPQ, breaks = c(0,10,20,30,40), 
                           labels = c("0-9", "10-19", "20-29", "30-40"), 
                           right = FALSE)
#split data frame by sampleID and then genome; summarise counts within each MAPQrange and return to single data frame
splitCount <- splitCount %>% 
  split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
      select(-sampleID) %>% 
      split(f = .$genome) %>% 
      lapply(function(z){z %>% 
          select(-genome, -MAPQ) %>% 
          group_by(MAPQrange) %>% 
          summarise_each(funs(sum))}) %>% 
      bind_rows(.id = "genome")}) %>% 
  bind_rows(.id = "sampleID")

#convert MAPQrange from factor to char?
splitCount <- splitCount %>% mutate_if(is.factor, as.character)

######### Does each genome have similar proportions of reads for each MAPQ?? ########

#Collate counts for rmdup for each MAPQrange
rmdupCount$MAPQrange <- cut(rmdupCount$MAPQ, 
                            breaks = c(0,10,20,30,40), 
                            labels = c("0-9", "10-19", "20-29", "30-40"), 
                            right = FALSE)
rmdupMAPQcount <- rmdupCount %>% 
  split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
      select(rmdupCount, MAPQrange) %>% 
      group_by(MAPQrange) %>% 
      summarise_each(funs(sum))}) %>% 
  bind_rows(.id = "sampleID")
rmdupMAPQcount <- rmdupMAPQcount %>% mutate_if(is.factor, as.character)

genome_names <- c(`aoris` = "A.oris", `cgracilis` = "C.gracilis", `esaphenum` = "E.saphenum", 
                  `fnucleatum` = "F.nucleatum", `mneoaurum` = "M.neourum", `moralis` = "M.oralis", 
                  `pgingivalis` = "P.gingivalis", `smutans` = "S.mutans", `tdenticola` = "T.denticola",
                  `tforsythia` = "T.forsythia")

#calculate proportion of total reads (for each MAPQrange) are represented by each genome & plot
pdf("plotsACADpres/initialData_MAPQbyGenome.pdf")
splitCount %>% 
  left_join(totalSplitCount, by = c("sampleID", "genome")) %>% 
  mutate(prop = splitCount.x/splitCount.y) %>% # calcualte prop for each genome at each MAPQrange
  split(f = .$sampleID) %>% #split df according to sampleID
  lapply(function(x){x %>% #use lapply to plot data for each sampleID separately
      ggplot(aes(x=MAPQrange, y=prop, fill=genome)) + 
      geom_bar(stat = "identity", position = position_dodge()) + 
      theme_bw() + 
      facet_wrap(~genome, ncol = 3, scales = "free", labeller = as_labeller(genome_names)) + 
      labs(x="MAPQ range", y="Proportion of total reads") + 
      scale_fill_manual(values = palette) + 
      guides(fill=FALSE) + 
      ggtitle(x$sampleID)})
dev.off()


############ Length Data #############

#create a list of lgdist files
lgDistFiles <- list.files("highQualAlnData", pattern = "lgdistribution.txt", 
                          full.names = TRUE, recursive = TRUE)
#read-in data from each text file and bind into a data frame
lengthData <- lgDistFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE, col_types = "cnn") %>%
    set_colnames(c("std", "length", "freq")) %>%
    mutate(fileName = x)
}) %>%
  bind_rows
#edit FileNames to include only the sampleID and genome
source("editFileNames.R")
lengthData$fileName <- editFileNames(lengthData)

lengthData$fileName <- gsub('highQualAlnData/results_', '', lengthData$fileName)
lengthData$fileName <- gsub('_Q30', '', lengthData$fileName)

#split fileName into sampleID and genome
lengthData <- colsplit(lengthData$fileName, "_", names = c("sampleID", "genome")) %>% 
  bind_cols(lengthData) %>% 
  select(-fileName)
#edit sampleID to be identical to previous
lengthData$sampleID <- gsub('ModernL7', 'Modern', lengthData$sampleID)
lengthData$sampleID <- gsub('ELSIDRON1', 'Elsidron1', lengthData$sampleID)

### Compare overall lengthData for each sample

#Collate lengths for each sample & plot distributions
#pdf("plotsACADpres/initialData_lengthDist_Q30.pdf")
lengthData %>% split(f = .$sampleID) %>%
  lapply(function(x){x %>% select(length, freq) %>% 
      group_by(length) %>% 
      summarise_each(funs(sum))}) %>% 
  bind_rows(.id = "sampleID") %>% 
    ggplot(aes(x=length, y=freq, colour=sampleID)) + 
    geom_line() + 
    theme_bw() + 
    labs(x="Read length", y="Number of reads", colour="Sample") + 
    ggtitle("Distribution of fragment lengths for each sample (Q>30)")
#dev.off()

#plot overall length data (boxplot)
expandedLengths <-  lengthData %>% 
  split(f = 1:nrow(.)) %>% 
  lapply(function(x){
    data_frame(sampleID = x$sampleID, 
               genome = x$genome, 
               length = rep(x$length, times = x$freq))}) %>% 
  bind_rows()

#BoxPlot
#pdf("plotsACADpres/initialData_lengthBoxPlot_Q30.pdf")
expandedLengths %>% ggplot(aes(x=sampleID, y=length, fill=sampleID)) + 
  geom_boxplot(outlier.color = "dark grey", outlier.size = 0.3) + 
  theme_bw() + theme(axis.title.x = element_blank()) + 
  ylab("Average Fragment Length") + guides(fill=FALSE) +
  ggtitle("Average fragment lengths per sample (Q > 30)")
#dev.off()

##Plot length distributions by genome
#pdf("plotsACADpres/initialData_lengthDistByGenome_Q30.pdf")
lengthData %>% 
  ggplot(aes(x=length, y=freq, colour=genome)) + 
  geom_line() + 
  theme_bw() + 
  labs(x="Read length", y="Number of reads", colour="Genome") + 
  facet_wrap(~sampleID) + 
  scale_colour_manual(values = palette,
                    name = "", 
                    labels=c(genome_names)) +
  ggtitle("Distribution of fragment lengths by genome (Q > 30)") 
#dev.off()

#Plot again, but facet_wrap genome
#pdf("plotsACADpres/initialData_lengthDistFacetByGenome_Q30.pdf")
lengthData %>% split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
    ggplot(aes(x=length, y=freq, colour=genome)) + 
    geom_line() + 
    theme_bw() + 
    labs(x="Read length", y="Number of reads", colour="Genome") + 
    facet_wrap(~genome, ncol=3, scales = "free_y", labeller = as_labeller(genome_names)) + 
    scale_colour_manual(values = palette) +
    theme(legend.position="none") + 
    ggtitle(x$sampleID, paste("(Q>30)"))})
#dev.off()

###Calculate key statistics about fragment lengths accroding to sample
lengthStatsBySampleQ30 <- expandedLengths %>% 
  select(sampleID, length) %>% 
  group_by(., sampleID) %>% 
  summarise(
    count = n(), 
    mean = mean(length, na.rm = TRUE), 
    sd = sd(length, na.rm = TRUE), 
    median = median(length, na.rm = TRUE), 
    IQR = IQR(length, na.rm = TRUE)
  )

################ phenotype characteristics #######################

phylum <- data_frame(genome = c("aoris", "cgracilis", "esaphenum", "fnucleatum",
                                "mneoaurum","moralis", "pgingivalis", "smutans", 
                                "tdenticola", "tforsythia"), 
                     phylum = c("Actinobacteria", "Proteobacteria", "Firmicutes", "Proteobacteria", 
                                "Actinobacteria", "Euryarchaeota", "Bacteroidetes", 
                                "Firmicutes", "Spirochaetes", "Bacteroidetes"))

GCcontent <- data_frame(genome = c("aoris", "cgracilis", "esaphenum", "fnucleatum",
                                   "mneoaurum","moralis", "pgingivalis", "smutans", 
                                   "tdenticola", "tforsythia"), 
                        GC = c("65-70", "45-50", "45-50", "<30", "65-70", "<30", "45-50",
                               "35-40", "35-40", "45-50"))

#Add information about bacterial cell wall
cellWall <- data_frame(genome = c("aoris", "cgracilis", "esaphenum", "fnucleatum", 
                                  "mneoaurum","moralis", "pgingivalis", "smutans", 
                                  "tdenticola", "tforsythia"), 
                       cellWall = c("Gram +", "Gram -", "Gram +", "Gram -", "Mycobacterium",
                                    "Archeal", "Gram -", "Gram +", "Gram -", "Gram -"))

######## Boxplots of Frag Length by Genome #########

#Bind this information to expandedLengths
expandedLengths <- expandedLengths %>% left_join(cellWall, by = "genome")
#Sort data according to cell wall
expandedLengths <- expandedLengths[order(expandedLengths$cellWall),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
expandedLengths$genome <- factor(expandedLengths$genome, levels = unique(expandedLengths$genome))

#boxplot Elsidron1
#pdf("plotsACADpres/initialData_ElsidronBoxPlot_Q30.pdf")
expandedLengths %>% filter(sampleID == "Elsidron1") %>% 
  ggplot(aes(x=genome, y=length, fill=cellWall)) + 
  geom_boxplot(outlier.colour = "dark grey", outlier.size = 0.3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(x="", y="Fragment length", fill="Cell wall") + 
  scale_x_discrete(labels=c(genome_names)) +
#  scale_y_log10() +
  scale_fill_manual(values=c(palette)) + 
  # annotate("text", x = c(1:10), y=5, size=3,
          # label = c("87,822", "37,037", "3,628", "7,115", "11,355", "29,099", "125,212", "8,356", "563", "1,011")) + 
  ggtitle("Elsidron 1 (Q>30)")
#dev.off()

#expandedLengths %>% filter(sampleID == "Elsidron1") %>% 
#  ggplot(aes(x = length, fill = cellWall)) +
#  geom_density() +
#  scale_x_log10() +
#  facet_wrap(~genome) +
#  theme_bw() +
#  theme(legend.position = c(0.75, 0.15))

#boxplot Modern
#pdf("plotsACADpres/initialData_ModernBoxPlot_Q30.pdf")
expandedLengths %>% filter(sampleID == "Modern") %>% 
ggplot(aes(x=genome, y=length, fill=cellWall)) + 
  geom_boxplot(outlier.colour = "dark grey", outlier.size = 0.3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(x="", y="Average fragment length", fill="Cell wall") + 
  scale_x_discrete(labels=c(genome_names)) +
  scale_fill_manual(values=c(palette)) + 
  annotate("text", x = c(1:10), y=5, size=3,
           label = c("185", "20,906", "272,789", "10,431", "6,760", "87,336", "17,540", "162", "1,285", "1,003")) + 
  ggtitle("Modern (Q>30)")
#dev.off()

####Add phylum and GC content data and repeat plots
#Bind this information to expandedLengths
expandedLengths <- expandedLengths %>% left_join(phylum, by = "genome")
#Sort data according to phylum
expandedLengths <- expandedLengths[order(expandedLengths$phylum),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
expandedLengths$genome <- factor(expandedLengths$genome, levels = unique(expandedLengths$genome))

expandedLengths %>% filter(sampleID == "Modern") %>% 
  ggplot(aes(x=genome, y=length, fill=phylum)) + 
  geom_boxplot(outlier.colour = "dark grey", outlier.size = 0.3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(x="", y="Average fragment length", fill="Phylum") + 
  scale_x_discrete(labels=c(genome_names)) +
  scale_fill_manual(values=c(palette)) + 
  ggtitle("Modern (Q>30)")

expandedLengths %>% filter(sampleID == "Elsidron1") %>% 
  ggplot(aes(x=genome, y=length, fill=phylum)) + 
  geom_boxplot(outlier.colour = "dark grey", outlier.size = 0.3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(x="", y="Average fragment length", fill="Phylum") + 
  scale_x_discrete(labels=c(genome_names)) +
  scale_fill_manual(values=c(palette)) + 
  ggtitle("Elsidron1 (Q>30)")

#Bind this information to expandedLengths
expandedLengths <- expandedLengths %>% left_join(GCcontent, by = "genome")
#Sort data according to GC
expandedLengths <- expandedLengths[order(expandedLengths$GC),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
expandedLengths$genome <- factor(expandedLengths$genome, levels = unique(expandedLengths$genome))

expandedLengths %>% filter(sampleID == "Elsidron1") %>% 
  ggplot(aes(x=genome, y=length, fill=GC)) + 
  geom_boxplot(outlier.colour = "dark grey", outlier.size = 0.3) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(x="", y="Average fragment length", fill="GC") + 
  scale_x_discrete(labels=c(genome_names)) +
  scale_fill_manual(values=c(palette)) + 
  ggtitle("Elsidron1 (Q>30)")

######### Repeated above using highQualAlnData - is there a difference??? ############

###Calculate key statistics about fragment lengths accroding to genome
lengthStatsByGenomeQ30 <- expandedLengths %>% 
  select(sampleID, length, genome) %>% 
  split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
  group_by(., genome) %>% 
  summarise(
    count = n(), 
    mean = mean(length, na.rm = TRUE), 
    sd = sd(length, na.rm = TRUE), 
    median = median(length, na.rm = TRUE), 
    IQR = IQR(length, na.rm = TRUE)
  )})

################## Substitution Rates ####################

#Create a vector with desired column names
subDataColNames <- (c("Chr", "End", "Std", "Pos", "A", "C", "G", "T", "Total", 
                      "GtoA", "CtoT", "AtoG", "TtoC", "AtoC", "AtoT", "CtoG", 
                      "CtoA", "TtoG", "TtoA", "GtoC", "GtoT"))

# create list of all .txt files in folder 
ntSubFiles <- list.files("highQualAlnData", pattern = "misincorporation.txt",
                         full.names = TRUE, recursive = TRUE)

# read-in files and bind into a data frame that include the FileName
ntSubData <- ntSubFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE, col_type = "cccnnnnnnnnnnnnnnnnnn---------") %>% 
    set_colnames(subDataColNames) %>% filter(Total > 0) %>% filter(Pos < 26) %>%
    mutate(FileName = x) %>% select(-Chr)
}) %>%
  bind_rows

#Edit FileName to include only the sampleID and Genome
ntSubData$FileName <- gsub('highQualAlnData/results_2NoAdapt_', '', ntSubData$FileName)
ntSubData$FileName <- gsub('_Q30_split/misincorporation.txt', '', ntSubData$FileName)
ntSubData$FileName <- gsub('ELSIDRON1', 'Elsidron1', ntSubData$FileName)

#Split FileName into separate columns - sampleID and Genome
ntSubData <- colsplit(ntSubData$FileName, "_", names = c("sampleID", "genome")) %>% 
  bind_cols(ntSubData) %>% 
  select(-FileName)

#Collate counts for pos and neg strand
ntSubData <- ntSubData %>% split(f = .$sampleID) %>%
  lapply(function(x){x %>% select(-sampleID) %>% split(f = .$End) %>% lapply(function(z){
  z %>% select(-End) %>% split(f = .$genome) %>% 
      lapply(function(a){a %>% select(-genome, -Std) %>% group_by(Pos) %>% summarise_each(funs(sum))}) %>% 
      bind_rows(.id = "genome")}) %>% 
      bind_rows(.id = "End")}) %>% bind_rows(.id = "sampleID")


################## Graphing ntSubData ##############


#Create an object containing genome counts for each sample
genomeCountsQ30 <- expandedLengths %>% 
  split(f = .$sampleID) %>% 
  lapply(function(x){x %>% 
      group_by(genome) %>% 
      summarise(count = n())}) %>% 
  bind_rows(.id = "sampleID")
#Convert Factor to character
genomeCountsQ30$genome <- as.character(genomeCountsQ30$genome)

#Bind this information to ntSubData
ntSubData <- ntSubData %>% left_join(genomeCountsQ30, by = c("sampleID", "genome"))

#Edit genome names (from aoris to A.oris)
source("editGenomeNames.R")
ntSubData$genome <- editGenomeNames(ntSubData)

# create graphing function
ntSub.graph <- function(df, na.rm = TRUE, ...){
  
  # Specify sampleID
  sampleID <- unique(df$sampleID)
  
  # create list of genomeID's in data to loop over 
  genomeID_list <- unique(df$genome)
  
  # create list of counts in data to loop over
  genomeCount_list <- unique(df$count)
  
  # create for loop to split data based on sampleID 
  for (i in seq_along(genomeID_list)) {
    
    # create object to store 5p data
    SubFreq_5p <- subset(df, df$genome==genomeID_list[i]) %>% filter(End == "5p") %>% 
      select(-End, -genome, -sampleID, -count) %>% 
      #group_by(Pos) %>% summarise_each(funs(sum)) %>% 
      mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
             AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
             TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
      select(-A, -C, -G, -T, -Total) %>% 
      melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
    
    #create object to store 3p data
    SubFreq_3p <- subset(df, df$genome==genomeID_list[i]) %>% filter(End == "3p") %>% 
      select(-End, -genome, -sampleID, -count) %>% 
      #group_by(Pos) %>% summarise_each(funs(sum)) %>% 
      mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
             AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
             TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
      select(-A, -C, -G, -T, -Total) %>% 
      melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
    
    #plot to object, 5p data
    plot5pData <- SubFreq_5p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + 
      geom_line() + 
      theme_bw() + 
      scale_y_continuous(limits = c(0, 0.55), position = "left") + 
      ylab("Substitution Frequency") + 
      xlab("Position from the 5' end") + 
      scale_colour_manual(values=c(palette), name = "Substitution") +
      theme(legend.position = "none")
    
    #plot to object, 3p data 
    plot3pData <- SubFreq_3p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + geom_line() + theme_bw() +
      theme(axis.title.y=element_blank()) + scale_x_reverse() + 
      scale_y_continuous(limits = c(0, 0.55), position = "right") +
      scale_colour_manual(values=c(palette), name = "Substitution") +
      ylab("Substitution Frequency") + xlab("Position from the 3' end")
    
    #print plots
    grid.arrange(plot5pData, plot3pData, ncol=2, widths=c(0.8,1), 
                 top = (paste(genomeID_list[i], '-',paste(sampleID,'(No. Reads = ', paste(genomeCount_list[i],')')))))
    
    #End loop
  }
}

# run graphing function on df
pdf("plotsACADpres/initialData_ElsidronNtSubQ30.pdf")
ntSubData %>% filter(sampleID == "Elsidron1") %>% ntSub.graph()
dev.off()

pdf("plotsACADpres/initialData_ModernNtSubQ30.pdf")
ntSubData %>% filter(sampleID == "Modern") %>% ntSub.graph()
dev.off()


############## comparison of C-T and G-A sub rate at pos 1 #####################

#Generate a list of 5pCtoT_freq.txt files
CtoT.Files <- list.files("highQualAlnData/", pattern = "5pCtoT_freq.txt", full.names = TRUE, recursive = TRUE)

#Create a function for reading in files and manipulating data; including adding phenotype characteristics
loadMisincorpData <- function(x){
  #load data
  data <- x %>% lapply(function(z){ 
    read_delim(z, delim = "\t", skip = 1, col_names = FALSE, col_types = cols("i", "n")) %>% # use lapply to read in each txt file
      set_colnames(c("pos", "freq")) %>% # specify the col_names
      mutate(fileName = z)}) %>% # add the fileName to the data
    bind_rows() # bind into single data frame
  #split fileName, retaining only sampleID and genome information in separate vectors
  data <- colsplit(data$fileName, "_", names = c("dir", "adapters", "sampleID", "genome", "extra")) %>% 
    bind_cols(data) %>% 
    select(-dir, -adapters, -fileName, -extra)
  #join cell wall information about each genome
  data <- data %>% left_join(cellWall, by = "genome")
  #join phylum information about each genome
  data <- data %>% left_join(phylum, by = "genome")
  #join GC content information about each genome
  data <- data %>% left_join(GCcontent, by = "genome")
  return(data)
  }

#Apply loadMisincorpData function to CtoT.Files
CtoTfreqData <- loadMisincorpData(CtoT.Files)
#Add information on genomeCounts
CtoTfreqData$sampleID <- gsub('ELSIDRON1', 'Elsidron1', CtoTfreqData$sampleID) 
CtoTfreqData <- genomeCountsQ30 %>% left_join(CtoTfreqData, by = c("sampleID", "genome"))

#Sort data based on sampleID [,1] then cellWall [,6]
CtoTfreqData <- arrange(CtoTfreqData, sampleID, cellWall)
#Convert genome to factor to prevent ggplot from re-ordering when plotting
CtoTfreqData$genome <- factor(CtoTfreqData$genome, levels = unique(CtoTfreqData$genome))

#Create function for plotting misincorp freq at Pos 1 to an object
plotMisFreq <- function(x, ...){x %>% ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size = 4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  ylab("Misincorporation frequency") + 
  xlab("") +
  ylim(0, 0.52) + 
  scale_colour_manual(values = palette) + 
  scale_x_discrete(labels = c(genome_names)) +  
  scale_shape_manual(values = c(19,17), name = "Sample ID")}

#Plot C->T freq at position 1 for all genomes and samples on the same axes
CtoTfreqData %>% 
  filter(pos == "1") %>% 
  plotMisFreq() + 
  scale_colour_manual(values = palette, name = "Cell wall") + 
  ggtitle("C -> T substitution frequency at position 1 5p' end")

#Repeat plot after removing genomes with < 500 reads
CtoTfreqData %>% filter(count > 500, pos == "1") %>% plotMisFreq() + labs(title = "C -> T substitution frequency at position 1 of 5' end",
                                                                               subtitle = "(No. Reads > 500)")
#Repeat plot after removing genomes with < 1000 reads
CtoTfreqData %>% filter(count > 1000, pos == "1") %>% plotMisFreq() + 
  scale_color_manual(values = palette, name = "Cell wall") + 
  labs(title = "C -> T misincorporation frequency at position 1 of 5' end",
      subtitle = "(No. Reads > 1000)")


### Read-in and plot G->A misincorporation
#Generate a list of 5pCtoT_freq.txt files
GtoA.Files <- list.files("highQualAlnData/", pattern = "3pGtoA_freq.txt", full.names = TRUE, recursive = TRUE)
#load and edit data from 3pGtoA_freq.text files using the editMisincorpData function
GtoAfreqData <- loadMisincorpData(GtoA.Files)
#Add information on readCount per genome (q > 30)
GtoAfreqData$sampleID <- gsub('ELSIDRON1', 'Elsidron1', GtoAfreqData$sampleID)
GtoAfreqData <- genomeCountsQ30 %>% left_join(GtoAfreqData, by = c("sampleID", "genome"))

#plot GtoA freq according to cellWall
GtoAfreqData <- arrange(GtoAfreqData, sampleID, cellWall)
GtoAfreqData$genome <- factor(GtoAfreqData$genome, levels = unique(GtoAfreqData$genome))
GtoAfreqData %>% filter(pos == "1") %>% plotMisFreq() + 
  scale_colour_manual(values = palette, name = "Cell wall") + 
  labs(title = "G -> A substitution frequency at position 1 of 3' end")

#plot GtoA freq according to cellWall, after filtering for reads > 1000
GtoAfreqData %>% filter(pos == "1", count > 1000) %>% plotMisFreq() + labs(title = "G -> A substitution frequence at position 1 of 3p' end", 
                                                                            subtitle = "(No. Reads > 1000)")


#plot GtoA freq according to phylum
GtoAfreqData <- arrange(GtoAfreqData, sampleID, phylum)
GtoAfreqData$genome <- factor(GtoAfreqData$genome, levels = unique(GtoAfreqData$genome))
GtoAfreqData %>% filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=phylum)) + geom_point(size=4) + 
  theme_bw() + 
  ylim(0, 0.52) + 
  scale_colour_manual(values = palette, name = "Phylum") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_x_discrete(labels = c(genome_names)) + 
  labs(y="Misincorporation frequency", x="", title="G -> A frequency, by phylum")

#plot CtoT freq according to phylum
CtoTfreqData <- arrange(CtoTfreqData, sampleID, phylum)
CtoTfreqData$genome <- factor(CtoTfreqData$genome, levels = unique(CtoTfreqData$genome))
CtoTfreqData %>% filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=phylum)) + geom_point(size=4) + 
  theme_bw() +
  ylim(0, 0.52) + 
  scale_x_discrete(labels = c(genome_names)) +
  scale_colour_manual(values = palette, name = "Phylum") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(y="Misincorporation frequency", x="", title="C -> T frequency, by phylum")
  

#plot GtoA freq according to GC content
GtoAfreqData <- arrange(GtoAfreqData, sampleID, GC)
GtoAfreqData$genome <- factor(GtoAfreqData$genome, levels = unique(GtoAfreqData$genome))
GtoAfreqData %>% filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=GC)) + 
  geom_point(size=4) + 
  theme_bw() + 
  scale_color_manual(values = palette, name = "GC content") + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(y="Misincorporation frequency", x="", title="C -> T frequency, by GC content")

#Plot CtoT freq according to GC content

################ Statistical Comparisons ###################

subFreq <- CtoTfreqData %>% filter(pos == "1") 
colnames(subFreq)[5] <- "CtoTfreq"
subFreq <- GtoAfreqData %>% filter(pos == "1") %>% left_join(subFreq)
colnames(subFreq)[5] <- "GtoAfreq"

### print length Data for comparison
head(expandedLengths)
#Add count data to expandedLengths
expandedLengths <- expandedLengths %>% left_join(genomeCountsQ30)

### print lengthStatsByGenome
lengthStatsByGenomeQ30

### print subFreq Data for comparison
subFreq

########## Questions ##############

# 1.Is there a statistically significant difference in the mean fragment lengths from each sample?

# 2. Is there a difference in the mean fragment lengths for each genome within a sample?

# 3. Is there a difference in mean fragment lengths for cellWall structure within a sample?

# 4.  Is there a difference in CtoT / GtoA sub frequency between samples?

# 5. Is there a difference in CtoT / GtoA sub frequency for cellWall structure within a sample?


##### Q1 kruskal.wallis - Null Hypothesis: There is no difference in mean fragment lengths between samples

##perform kruskal.wallis test on expandedLengths
#must first convert sampleID from chr to factor
expandedLengths$sampleID <- factor(expandedLengths$sampleID)
kruskal.test(length ~ sampleID, data = expandedLengths) 
  #Indicates there is a statistically significant difference


##### Q2 kruskal.wallis - Null Hypothesis: There is no difference in mean fragment lengths between genomes

#Split the expandedLengths data frame according to sampleID into separate data frames, each with name of the sample
list2env(split(expandedLengths, expandedLengths$sampleID), envir = .GlobalEnv)

#Convert genome from chr to factor
Modern$genome <- factor(Modern$genome, levels = unique(Modern$genome))

##perform kruskal.wallis test on Modern Lengths
kw.ModernLengths <- kruskal.test(length ~ genome, data = Modern)
#Indicates there is a statistically significant difference in lengths.

##remove genomes with low number of reads aligning and re-run kruskal.test
Modern.count1000 <- Modern %>% filter(count > 1000)
##re-run kruskal.wallis test
kw.ModernLengths.count1000 <- kruskal.test(length ~ genome, data = Modern.count1000) # still indicates a statistically sig difference

###Dunn's test?? / Conover-Iman test ?? (to identify genomes which demonstrate significantly different lengths)



#### Q1 perm.test - Null Hypothesis: There is no difference in mean fragment lengths of each sample

###Permutation tests
library(coin)
str(expandedLengths)
##Permutation test of independence
#This test treats the two groups (e.g. Elsidron1 & Modern) as independent samples, 
  #and tests if there is a difference in values between the two groups. 
independence_test(length ~ sampleID, data = expandedLengths)

##Permutation test of symmetry
#This test treats the two groups (e.g. Elsidron1 & Modern) as having paired or repeated data, paired within sample. 

#### Q4 perm.test - Null Hypothesis: There is no difference in the C->T frequency observed between samples
#remove GtoAfreq data
subFreqCtoT <- subFreq %>% select(-GtoAfreq)
#convert sampleID from chr to factor
subFreqCtoT$sampleID <- factor(subFreqCtoT$sampleID)
str(subFreqCtoT)
#perform independence_test
independence_test(CtoTfreq ~ sampleID, data = subFreqCtoT) #indicates significant difference, but less sig than between frag lengths??

#### Q2 perm.test - Null Hypothesis: There is no difference in mean fragment lengths between genomes
str(Modern.count1000)
independence_test(length ~ genome, data = Modern.count1000)

################# lmem done with data from 6 samples of 15 genomes ###################

##Collate genome counts for Expanded CtoT data
GenomeCountExpandedQ30 <- read_delim("mapData/highQualMapData/high_qual_split_count.txt", delim = "\t", 
           skip = 1, col_names = FALSE, col_types = "cn") %>%
  set_colnames(c("fileName", "genomeCount"))
#Edit fileName to include only sampleID and genome
GenomeCountExpandedQ30$fileName <- editFileNames(GenomeCountExpandedQ30)
#Split sampleID and Genome into separate vector
GenomeCountExpandedQ30 <- colsplit(GenomeCountExpandedQ30$fileName, "_", names = c("sampleID", "genome")) %>% 
  bind_cols(GenomeCountExpandedQ30) %>% select(-fileName)

#Generate vectors with information about cell wall type, GC content, oxidase tests, phylum
cellWallExpanded <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                  "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                  "T.denticola", "T.forsythia"), 
                       cellWall = c("GramPos", "GramPos", "GramNeg", "GramPos", "GramNeg", "GramNeg", "Mycobacterium",
                                    "Archeal", "GramNeg", "GramNeg", "GramNeg", "GramPos", "GramPos", "GramNeg", "GramNeg"))

GCcontentExpanded <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                   "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                   "T.denticola", "T.forsythia"), 
                        GC = c("65-70", "45-50", "45-50", "45-50", "<30", "35-40", "65-70", "<30", "50-55", "45-50", "40-45",
                               "40-45", "35-40", "35-40", "45-50"))

phylumExpanded <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                "T.denticola", "T.forsythia"), 
                     phylum = c("Actinobacteria", "Actinobacteria", "Proteobacteria", "Firmicutes", "Fusobacteria", "Proteobacteria", 
                                "Actinobacteria", "Euryarchaeota", "Proteobacteria", "Bacteroidetes", "Bacteroidetes", "Firmicutes", 
                                "Firmicutes", "Spirochaetes", "Bacteroidetes"))

catalaseExpanded <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                  "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                  "T.denticola", "T.forsythia"), 
                       oxidase = c("V", "N", "N", "N", "N", "P", "P", "U", "P", "N", "N", "N", "N", "N", "V"))

#Generate a list of 5pCtoT_freq.txt files
CtoT.Files.Expanded <- list.files("mapData/highQualMapData", pattern = "5pCtoT_freq.txt", full.names = TRUE, recursive = TRUE)
# exclude SpyNEW as too few reads to analyse damage patterns
CtoT.Files.Expanded <- CtoT.Files.Expanded[!grepl("SPYNEWL8", CtoT.Files.Expanded)]
#Read in each of these files and bind together in a data frame (MUST specify col_types!)
CtoTfreqDataExpanded <- CtoT.Files.Expanded %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 1, col_names = FALSE, col_types = cols("i", "n")) %>% 
    set_colnames(c("pos", "freq")) %>% 
    mutate(fileName = x)}) %>% 
  bind_rows()
#Edit fileName to include only sampleID and genome
CtoTfreqDataExpanded$fileName <- editFileNames(CtoTfreqDataExpanded)
#Split sampleID and Genome into separate vector
colsplit(CtoTfreqDataExpanded$fileName, "_", names = c("sampleID", "genome")) %>% 
  bind_cols(CtoTfreqDataExpanded) %>% 
  select(-fileName) -> CtoTfreqDataExpanded
#Add information on cell.wall and bind to data frame
CtoTfreqDataExpanded <- CtoTfreqDataExpanded %>% left_join(cellWallExpanded, by = "genome")
#Sort data based on sampleID [,1] then cellWall [,5]
CtoTfreqDataExpanded <- CtoTfreqDataExpanded[order( CtoTfreqDataExpanded[,1], CtoTfreqDataExpanded[,5] ), ]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
CtoTfreqDataExpanded$genome <- factor(CtoTfreqDataExpanded$genome, levels = unique(CtoTfreqDataExpanded$genome))

#Plot all data on same axes
CtoTfreqDataExpanded %>% 
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ylim(0, 0.52) + 
  #  theme(legend.position = "none") +
  ggtitle("Cytosine to Thymine misincorporation frequency\nat Position 1 of the 5' end")

###Add information on number of reads aligning to each sample (from genomeCount) and filter out genomes with <1000 reads
CtoTfreqDataExpanded <- left_join(CtoTfreqDataExpanded, GenomeCountExpandedQ30, by = c("sampleID", "genome"))
#Sort data based on sampleID [,1] then cellWall [,5]
CtoTfreqDataExpanded <- CtoTfreqDataExpanded[order( CtoTfreqDataExpanded[,1], CtoTfreqDataExpanded[,5] ), ]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
CtoTfreqDataExpanded$genome <- factor(CtoTfreqDataExpanded$genome, levels = unique(CtoTfreqDataExpanded$genome))
#Filter out genomes with < 1000 reads and plot
CtoTfreqDataExpanded %>% 
  filter(genomeCount > 1000) %>%
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ylim(0, 0.52) + 
  #  theme(legend.position = "none") +
  ggtitle("Cytosine to Thymine misincorporation frequency\nat Position 1 of the 5' end (Reads > 1000)")  

#Generate object containing CtoT freq only for pos 1, with any 0 values removed
pos1CtoTallSamples <- CtoTfreqDataExpanded %>% filter(pos == "1", freq != "0")

pos1CtoTallSamples %>% 
  filter(!sampleID %in% c("Modern", "CHIMP")) %>% 
  filter(!cellWall %in% c("Archeal", "Mycobacterium")) %>% 
  ggplot(aes(x=genome, y=freq, colour=sampleID, shape=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ylim(0, 0.52) + 
  ggtitle("C -> T frequency at position 1, excluding 0 values")

#perform linear mixed model effects on CtoTfreq by cellWall, nested within ancient samples
lmer(freq~cellWall + (1|sampleID), data = filter(pos1CtoTallSamples, grepl("Gram", cellWall), 
                                                 !sampleID %in% c("Modern", "CHIMP"))) %>% summary()

abunFilteredpos1CtoTallSamples %>% 
  filter(!sampleID %in% c("Modern", "CHIMP", "SPYOLD")) %>% 
  filter(!cellWall %in% c("Archeal", "Mycobacterium")) %>% 
  ggplot(aes(x=genome, y=freq, colour=sampleID, shape=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ylim(0, 0.52) + 
  ggtitle("C -> T frequency at position 1, excluding low abundand genomes")  

#Filter out low abundant genomes and repeat lmer
abunFilteredpos1CtoTallSamples <- pos1CtoTallSamples %>% filter(genomeCount > 1000)
lmer(freq~cellWall + (1|sampleID), data = filter(abunFilteredpos1CtoTallSamples, grepl("Gram", cellWall), 
                                                 !sampleID %in% c("Modern", "CHIMP"))) %>% summary()
#Again without SpyOld
lmer(freq~cellWall + (1|sampleID), data = filter(abunFilteredpos1CtoTallSamples, grepl("Gram", cellWall), 
                                                 !sampleID %in% c("Modern", "CHIMP", "SPYOLD"))) %>% summary()


#density plot of log transformed lengths for Elsidron1
expandedLengths %>% filter(sampleID == "Elsidron1") %>% 
   ggplot(aes(x = length, fill = cellWall)) +
   geom_density() +
   scale_x_log10() +
   facet_wrap(~genome) +
   theme_bw() +
   theme(legend.position = c(0.75, 0.15))

#calculate mean_log_lengths and assign to new object
expandedLengths %>%
   group_by(sampleID, genome) %>%
   summarise(mn_length = mean(logLength)) %>%
   rename(name = genome) %>%
   mutate(genome = str_replace_all(name, "\\.", ""),
          genome = tolower(genome)) %>%
   left_join(cellWall) %>%
   mutate(Modern = sampleID %in% c("Modern")) %>%
   filter(!is.na(cellWall),
                  grepl("Gram", cellWall), sampleID != "CHIMP") -> tempLengths

#perform linear mixed model effects on mn_log_lengths by cell wall, nested in ancient samples
summary(lmer(mn_length~cellWall + (1|sampleID), data = tempLengths, 
             subset = !sampleID %in% c("Modern", "CHIMP")))

#plot of CtoTfreq against GC (with cellWall as colour)
CtoTfreqData %>% ggplot(aes(x = GC, y = freq, colour = cellWall)) + 
  geom_point() + 
  facet_wrap(~sampleID, scales = "free_y") + 
  theme(legend.position = c(0.8, 0.25))


#perform linear mixed model effects on CtoTfreq by cellWall, nested within ancient samples
lmer(freq~cellWall + (1|sampleID), data = filter(CtoTfreqData, grepl("Gram", cellWall), 
                                                 !sampleID %in% c("Modern", "CHIMP"))) %>% summary()

#as above but with 0 frequencies removed (***Must check why 0 exist in data!!! error during import??)
lmer(freq~cellWall + (1|sampleID), data = filter(CtoTfreqData, grepl("Gram", cellWall), 
                                                 !sampleID %in% c("Modern", "CHIMP"), freq != 0)) %>% 
  summary()
#convert t value to p-value using df = 7
library(lmerTest)
pt(-2.5, df = 7)

#attempts to include both sampleID and cellWall as contributing to variance
lm(freq~(sampleID + cellWall)^2 , data = filter(CtoTfreqData, 
                                                !sampleID %in% c("Modern", "CHIMP")), 
   subset = grepl("Gram", cellWall)) %>% summary
#only some samples (e.g. SPYOLD) indicate significant difference (as affected by contamination)
#cell wall is different

#model using sampleID and cellWall as effects + GC
lm(freq~(sampleID + cellWall)^2 + GC , data = filter(CtoTfreqData, 
                                                     !sampleID %in% c("Modern", "CHIMP")), 
   subset = grepl("Gram", cellWall)) %>% summary
#GC does not appear to be having an effect.
#SPYOLD sample significantly different to others (p-value =3.67e-08)
#cellWall p-value=0.0215

#Model using sample, cell wall and GC as effects
lm(freq~(sampleID + cellWall + GC)^2 , data = filter(CtoTfreqData, 
                                                     !sampleID %in% c("Modern", "CHIMP")), 
   subset = grepl("Gram", cellWall)) %>% summary
#only SPYOLD and cell wall sig different


#attempts to see if GC content inlfuencing CtoT freq


###################### Running on single genome - S.mutans #######################

singleGenomeReadCounts <- read.csv(file="smutans/read_count.txt", sep="", 
                                    skip = 1, header = FALSE, col.names = c("count", "MAPQ"))
singleGenomeReadCounts <- singleGenomeReadCounts %>% 
  mutate(bam = grepl("bam", count), fileNo = cumsum(bam)) %>% split(f = .$fileNo)
names(singleGenomeReadCounts) <- c("bwa", "rmdup", "sorted")
singleGenomeReadCounts <- singleGenomeReadCounts %>% 
  bind_rows(.id = "fileName") %>% 
  select(-bam, -fileNo) %>% filter(MAPQ != "NA", fileName == "rmdup") %>%
  mutate(genome = "smutans")

##Compare counts to those in 10 alignment and in 15 alignment
smutansSplitCount <- splitCount %>% filter(sampleID == "ELSIDRON1L7", genome == "smutans") %>% select(-sampleID)
#combint both counts
smutansSplitCount %>% left_join(singleGenomeReadCounts, by = c("genome", "MAPQ"))

##Find data on 15 genome counts
splitCount_15genomes <- read.csv(file = "mapData/lowQualMapData/low_qual_split_count.txt", sep = "", 
                                 skip = 1, header = FALSE, col.names = c("count", "MAPQ"))
splitCount_15genomes <- splitCount_15genomes %>% 
  mutate(bam = grepl("bam", count), fileNo = cumsum(bam)) %>% split(f = .$fileNo)
splitCount_15_Files <- list.files(path = "mapData/lowQualMapData", pattern = "_split.bam")
names(splitCount_15genomes) <- splitCount_15_Files
splitCount_15genomes <- splitCount_15genomes %>% bind_rows(.id = "fileName") %>% select(-bam, -fileNo) %>% filter(MAPQ != "NA")
#edit fileName
splitCount_15genomes$fileName <- editFileNames(splitCount_15genomes)
splitCount_15genomes$fileName <- gsub('_split.bam', '', splitCount_15genomes$fileName)
smutansSplitCount_15Genomes <- splitCount_15genomes %>% filter(fileName == "ELSIDRON1_S.mutans")
smutansSplitCount_15Genomes <- smutansSplitCount_15Genomes %>% mutate(genome = "smutans") %>% select(-fileName)
names(smutansSplitCount_15Genomes) <- c("comp15Count", "MAPQ", "genome")
smutansCounts <- smutansSplitCount %>% left_join(singleGenomeReadCounts) %>% left_join(smutansSplitCount_15Genomes)
smutansCounts$count <- as.numeric(as.character(smutansCounts$count))
smutansCounts$comp15Count <- as.numeric(as.character(smutansCounts$comp15Count))

#Add categorical MAPQranges and plot just 10 and 15 genomes
smutansCounts <- smutansCounts %>% select(-genome, -fileName)
smutansCounts$MAPQrange <- cut(smutansCounts$MAPQ, breaks = c(0,10,20,30,40), 
                           labels = c("0-9", "10-19", "20-29", "30-40"), 
                           right = FALSE)
smutansCounts <- smutansCounts %>% select(-MAPQ) %>% group_by(MAPQrange) %>% summarise_each(funs(sum))

smutansCounts %>% select(-count) %>% 
  melt(id.vars = "MAPQrange") %>% 
  ggplot(aes(x=MAPQrange, y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  scale_fill_manual(values = palette, name = "Alignment", labels = c("No other Strep. species", "Alignment with S. mitis")) + 
  labs(x="MAPQ range", y="Number of alignments", title="Effect of alignment against related species (S. mutans)")







