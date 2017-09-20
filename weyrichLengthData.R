#Comparison of length data for genomes aligning to Elsidron1 and Modern

#load packages
library(dplyr)
library(readr)
library(magrittr)
library(tibble)
library(stringr)
library(reshape2)
library(ggplot2)

# create list of all .txt files in folder 
lgDistFiles <- list.files("alnData", pattern = "lgdistribution.txt", 
                          full.names = TRUE, recursive = TRUE)
lengthData <- lgDistFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>%
    set_colnames(c("Std", "Length", "Occurences")) %>%
    mutate(FileName = x)
}) %>%
  bind_rows

#edit FileNames to include only the sampleID and genome
lengthData$FileName <- gsub('alnData/results_2NoAdapt_', '', lengthData$FileName)
lengthData$FileName <- gsub('_split/lgdistribution.txt', '', lengthData$FileName)
lengthData$FileName <- gsub('L7_lTACTG_rCTCGA', '', lengthData$FileName)
lengthData$FileName <- gsub('L7_lAAGAG_rNONE', '', lengthData$FileName)
lengthData$FileName <- gsub('alnData/results_', '', lengthData$FileName)
lengthData$FileName <- gsub('L7', '', lengthData$FileName)
lengthData$FileName <- gsub('ELSIDRON1', 'Elsidron1', lengthData$FileName)

#Reshape data to separate sampleID and genome into separate vectors
lengthData <- lengthData %>% 
  mutate(Genome = str_extract(FileName, "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)"), 
         sampleID = str_replace(FileName, "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)", ""))

#sort columns so that genomes listed in alphabetical order
lengthData <- lengthData[with(lengthData, order(sampleID, Genome)), ]

# create graphing function
lgDist.graph <- function(df, na.rm = TRUE, ...){
  
  # create list of sampleID's in data to loop over 
  sampleID_list <- unique(df$sampleID)
  
  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(sampleID_list)) {

    # create histogram plot for each sample in df 
    plot1 <-
      ggplot(subset(df, df$sampleID==sampleID_list[i]),
        aes(Length, Occurences, fill = Std)) + 
      geom_bar(stat = "identity") +
      facet_wrap(( ~ Genome), ncol=4) +
      theme_bw() +
      theme(legend.position="none") + 
      scale_y_continuous("Frequency") +
      scale_x_continuous("Fragment length") +
      ggtitle(paste(sampleID_list[i], ' Fragment Lengths by Genome'))
    
    plot2 <-
      ggplot(subset(df, df$sampleID==sampleID_list[i]),
             aes(Length, Occurences, colour = Genome)) +
      geom_line() +
      xlab("Fragment length") +
      ggtitle(paste(sampleID_list[i], ' Fragment Lenght Distribution')) +
      theme_bw() 
    
    print(plot1)
    print(plot2)
    
  }
} 

# run graphing function on long df
lgDist.graph(lengthData) 


################ Compare average fragment length of genomes #######################

#use lapply to read each text file, convert to a vector of lengths and bind together
lengthData2 <- lgDistFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>%
    select(-1) %>% group_by(X2) %>% summarise_each(funs(sum)) -> length.vector
  rep(length.vector$X2, length.vector$X3)
}) 
#Add names to each object in the list
names(lengthData2) <- lgDistFiles
#Bind list into a dataframe with FileName listed
lengthVectors_long <- names(lengthData2) %>% 
  lapply(function(x){data_frame(FileName = x, lengths = lengthData2[[x]])}) %>% bind_rows()

#edit FileNames to include only the sampleID and genome
lengthVectors_long$FileName <- gsub('alnData/results_2NoAdapt_', '', lengthVectors_long$FileName)
lengthVectors_long$FileName <- gsub('_split/lgdistribution.txt', '', lengthVectors_long$FileName)
lengthVectors_long$FileName <- gsub('L7_lTACTG_rCTCGA', '', lengthVectors_long$FileName)
lengthVectors_long$FileName <- gsub('L7_lAAGAG_rNONE', '', lengthVectors_long$FileName)
lengthVectors_long$FileName <- gsub('alnData/results_', '', lengthVectors_long$FileName)
lengthVectors_long$FileName <- gsub('L7', '', lengthVectors_long$FileName)
lengthVectors_long$FileName <- gsub('ELSIDRON1', 'Elsidron1', lengthVectors_long$FileName)

#Reshape data to separate sampleID and genome into separate vectors
lengthVectors_long <- lengthVectors_long %>%
  mutate(Genome = str_extract(FileName, "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)"), 
         sampleID = str_replace(FileName, "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)", ""))

#sort columns so that genomes listed in alphabetical order
lengthVectors_long <- lengthVectors_long[with(lengthVectors_long, order(sampleID, Genome)), ]
#remove FileName vector from data frame
#lengthVectors_long <- lengthVectors_long %>% select(-FileName)
#Add information about bacterial cell wall
cell.wall <- data_frame(Genome = c("aoris", "cgracilis", "esaphenum", "fnucleatum", 
                                   "mneoaurum","moralis", "pgingivalis", "smutans", 
                                   "tdenticola", "tforsythia"), 
                       cell.wall = c("gram+", "gram-", "gram+", "gram-", "Mycobacterium",
                                     "Archeal", "gram-", "gram+", "gram-", "gram-"))

#join this information to overall dataframe with left_join
lengthVectors_long <- lengthVectors_long %>% left_join(cell.wall, by = "Genome")

# create graphing function
length.boxplot <- function(df, na.rm = TRUE, ...){
  
  # create list of sampleID's in data to loop over 
  sampleID_list <- unique(df$sampleID)
  
  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(sampleID_list)) {
    
    # create histogram plot for each sample in df 
    boxplot1 <-
      ggplot(subset(df, df$sampleID==sampleID_list[i]),
             aes(x=Genome, y=lengths, fill=cell.wall)) + 
      stat_boxplot(geom = 'errorbar') + geom_boxplot(outlier.shape = NA) +
      theme_bw() + ylim(30, 125) +
      ylab("Average fragment length") + xlab("") +
      scale_x_discrete(labels=c("A.oris", "C.gracilis", "E.saphenum", "F.nucleatum", 
              "M.neoaurum", "M.oralis", "P.gingivalis", "S.mutans", "T.denticola", 
              "T.forsythia")) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
      ggtitle(paste(sampleID_list[i], ' Sample', "\nAverage fragment lengths by genome"))
    
    print(boxplot1)
  }
}

# run graphing function on long df
length.boxplot(lengthVectors_long)

#Produce summary statistics on Elsidron1 data
split(lengthVectors_Elsidron$lengths, lengthVectors_Elsidron$Genome) %>% lapply(function(x){summary(x)})

#Generate a function to apply to stat_summary so that can plot max and min values
# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# do it - generate plot
lengthVectors_long %>% 
  ggplot(aes(x=Genome, y=lengths, fill=cell.wall)) + stat_summary(fun.data = f, geom="boxplot") +
  theme_bw() + ylab("Average fragment length") + xlab("") +
  scale_x_discrete(labels=c("A.oris", "C.gracilis", "E.saphenum", "F.nucleatum", 
                            "M.neoaurum", "M.oralis", "P.gingivalis", "S.mutans", "T.denticola", 
                            "T.forsythia")) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) + facet_wrap(~sampleID, ncol=1)

########################### Length Plots using Q30 data ############################

# create list of all .txt files in folder 
lgDistFilesQ30 <- list.files("highQualAlnData", pattern = "lgdistribution.txt", 
                          full.names = TRUE, recursive = TRUE)
# Read-in data from each file and bind into a single data frame
lengthDataQ30 <- lgDistFilesQ30 %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>%
    set_colnames(c("Std", "Length", "Occurences")) %>%
    mutate(FileName = x)
}) %>%
  bind_rows

#edit FileNames to include only the sampleID and genome
lengthDataQ30$FileName <- gsub('highQualAlnData/results_2NoAdapt_', '', lengthDataQ30$FileName)
lengthDataQ30$FileName <- gsub('_Q30_split/lgdistribution.txt', '', lengthDataQ30$FileName)
lengthDataQ30$FileName <- gsub('ELSIDRON1', 'Elsidron1', lengthDataQ30$FileName)

#Reshape data to separate sampleID and genome into separate vectors
lengthDataQ30 <- lengthDataQ30 %>% 
  mutate(Genome = str_extract(FileName, 
        "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)"), 
         sampleID = str_replace(FileName, 
        "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)", ""))

# create graphing function
Q30.lgDist.graph <- function(df, na.rm = TRUE, ...){
  
  # create list of sampleID's in data to loop over 
  sampleID_list2 <- unique(df$sampleID)
  
  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(sampleID_list2)) {
    
    # create histogram plot for each sample in df 
    Q30lgDist.plot1 <-
      ggplot(subset(df, df$sampleID==sampleID_list2[i]),
             aes(Length, Occurences, fill = Std)) + 
      geom_bar(stat = "identity") +
      facet_wrap(( ~ Genome), ncol=4) +
      theme_bw() +
      theme(legend.position="none") + 
      scale_y_continuous("Frequency") +
      scale_x_continuous("Fragment length") +
      ggtitle(paste(sampleID_list2[i], ' Fragment Lengths by Genome'))
    
    Q30lgDist.plot2 <-
      ggplot(subset(df, df$sampleID==sampleID_list2[i]),
             aes(Length, Occurences, colour = Genome)) +
      geom_line() +
      xlab("Fragment length") +
      ggtitle(paste(sampleID_list2[i], ' Fragment Lenght Distribution')) +
      theme_bw() 
    
    print(Q30lgDist.plot1)
    print(Q30lgDist.plot2)
    
  }
} 

# run graphing function on long df
Q30.lgDist.graph(lengthDataQ30) 

### Compare average fragment length of genomes

#use lapply to read each text file, convert to a vector of lengths and bind together
Q30lengthVectors <- lgDistFilesQ30 %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>%
    select(-1) %>% group_by(X2) %>% summarise_each(funs(sum)) -> length.vector
  rep(length.vector$X2, length.vector$X3)
}) 
#Add names to each object in the list
names(Q30lengthVectors) <- lgDistFilesQ30
#Bind list into a dataframe with FileName listed
Q30lengthVectors_long <- names(Q30lengthVectors) %>% 
  lapply(function(x){data_frame(FileName = x, lengths = Q30lengthVectors[[x]])}) %>% bind_rows()

#edit FileNames to include only the sampleID and genome
Q30lengthVectors_long$FileName <- gsub('highQualAlnData/results_2NoAdapt_', '', Q30lengthVectors_long$FileName)
Q30lengthVectors_long$FileName <- gsub('_Q30_split/lgdistribution.txt', '', Q30lengthVectors_long$FileName)

#Reshape data to separate sampleID and genome into separate vectors
Q30lengthVectors_long <- Q30lengthVectors_long %>%
mutate(Genome = str_extract(FileName, 
                  "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)"), 
       sampleID = str_replace(FileName, 
             "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)", ""))
#remove FileName vector from data frame
Q30lengthVectors_long <- Q30lengthVectors_long %>% select(-FileName)
#Add information about bacterial cell wall
cell.wall <- data_frame(Genome = c("aoris", "cgracilis", "esaphenum", "fnucleatum", 
                                   "mneoaurum","moralis", "pgingivalis", "smutans", 
                                   "tdenticola", "tforsythia"), 
                        cell.wall = c("gram+", "gram-", "gram+", "gram-", "Mycobacterium",
                                      "Archeal", "gram-", "gram+", "gram-", "gram-"))

#join this information to overall dataframe with left_join
Q30lengthVectors_long <- Q30lengthVectors_long %>% left_join(cell.wall, by = "Genome")

# create histogram plot for each sample
Q30.boxplot1 <- Q30lengthVectors_long %>% ggplot(aes(x=Genome, y=lengths, fill=cell.wall)) + 
  stat_boxplot(geom = 'errorbar') + geom_boxplot(outlier.shape = NA) +
  theme_bw() + ylim(30, 125) +
  ylab("Average fragment length") + xlab("") +
  scale_x_discrete(labels=c("A.oris", "C.gracilis", "E.saphenum", "F.nucleatum", 
                            "M.neoaurum", "M.oralis", "P.gingivalis", "S.mutans", "T.denticola", 
                            "T.forsythia")) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  ggtitle("Average fragment lengths by genome") + facet_wrap(~sampleID)



################## Normalize counts based on sample size ##########################

#Read in split_aligned_read_count.txt and assign to object splitCount
splitCount <- read_delim(file = "alnData/split_aligned_read_count.txt", delim = "\t", skip = 1, col_names = FALSE) %>% 
  set_colnames(c("FileName", "totReadCount"))
#Edit FileName to include only sampleID and genome
splitCount$FileName <- gsub('2NoAdapt_', '', splitCount$FileName)
splitCount$FileName <- gsub('_split/lgdistribution.txt', '', splitCount$FileName)
splitCount$FileName <- gsub('L7_lTACTG_rCTCGA', '', splitCount$FileName)
splitCount$FileName <- gsub('_lAAGAG_rNONE', '', splitCount$FileName)
#Edit FileName to ensure sampleID is consistent
splitCount$FileName <- gsub('Elsidron1', 'ELSIDRON1', splitCount$FileName)
#Join splitCount data with lengthData to create normLengthData
normLengthData <- lengthData %>% left_join(splitCount, by = "FileName")
#normalize occurences of each read length by converting to proportion of sample
normLengthData <- normLengthData %>% mutate(Prop = (Occurences/totReadCount)*100)

##### Plot normalised Length distribution for each genome in a sample ####

#Plot for Elsidron1 sample
normLengthData %>% filter(sampleID == "ELSIDRON1") %>% 
  ggplot(aes(x=Length, y=Prop, colour = Genome)) + geom_line() + theme_bw() +
  xlab("Fragment length") + ylab("Proportion of reads") + ggtitle("Fragment length distribution for Elsidron1")
#Plot for ModernL7 sample
normLengthData %>% filter(sampleID == "ModernL7") %>% 
  ggplot(aes(x=Length, y=Prop, colour = Genome)) + geom_line() + theme_bw() +
  xlab("Fragment length") + ylab("Proportion of reads") + ggtitle("Fragment length distribution for ModernL7")
#Plot for ModernL7, excluding M.oralis
normLengthData %>% filter(sampleID == "ModernL7") %>% filter(Genome != "moralis") %>% 
  ggplot(aes(x=Length, y=Prop, colour = Genome)) + geom_line() + theme_bw() +
  xlab("Fragment length") + ylab("Proportion of reads") + 
  ggtitle("Fragment length distribution for ModernL7, excluding M.oralis")

#### plot normalised histogram data for each sample, and genome

#Plot for Elsidron1
normLengthData %>% filter(sampleID == "ELSIDRON1") %>% 
  ggplot(aes(x=Length, y=Prop, fill=Std)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_wrap((~Genome), ncol = 4) + theme(legend.position="none") + 
  scale_y_continuous("Proportion of total reads") +
  scale_x_continuous("Fragment length") +
  ggtitle("Fragment Lengths by Genome for Elsidron1")
  
#Plot for ModernL7
normLengthData %>% filter(sampleID != "ELSIDRON1") %>% 
  ggplot(aes(x=Length, y=Prop, fill=Std)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_wrap((~Genome), ncol = 4) + theme(legend.position="none") + 
  scale_y_continuous("Proportion of total reads") +
  scale_x_continuous("Fragment length") +
  ggtitle("Fragment Lengths by Genome for ModernL7")

#Plot for ModernL7, excluding M.oralis
normLengthData %>% filter(sampleID != "ELSIDRON1") %>% filter(Genome != "moralis") %>%
  ggplot(aes(x=Length, y=Prop, fill=Std)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_wrap((~Genome), ncol = 4) + theme(legend.position="none") + 
  scale_y_continuous("Proportion of total reads") +
  scale_x_continuous("Fragment length") +
  ggtitle("Fragment Lengths by Genome for ModernL7, excluding M.oralis")

############# Plot normalized Q30 fragment lengths #####################

#Read in highQual_split_aligned_read_count.txt and assign to object splitCount
highQualSplitCount <- read_delim(file = "highQualAlnData/highQual_split_aligned_read_count.txt", delim = "\t",
                                 skip = 1, col_names = FALSE) %>% set_colnames(c("FileName", "totReadCount"))

#Edit FileName to include only sampleID and genome
highQualSplitCount$FileName <- gsub('2NoAdapt_', '', highQualSplitCount$FileName)

#Edit FileName to ensure sampleID is consistent
highQualSplitCount$FileName <- gsub('ELSIDRON1', 'Elsidron1', highQualSplitCount$FileName)

#Join splitCount data with lengthData to create normLengthData
normLengthDataQ30 <- lengthDataQ30 %>% left_join(highQualSplitCount, by = "FileName")

#normalize occurences of each read length by converting to proportion of sample
normLengthDataQ30 <- normLengthDataQ30 %>% mutate(Prop = (Occurences/totReadCount)*100)

###Plot disribution data

#Plot for Elsidron1 sample
normLengthDataQ30 %>% filter(sampleID == "Elsidron1") %>% 
  ggplot(aes(x=Length, y=Prop, colour = Genome)) + geom_line() + theme_bw() +
  xlab("Fragment length") + ylab("Proportion of reads") + ggtitle("Fragment length distribution for Elsidron1")
#Plot for ModernL7 sample
normLengthDataQ30 %>% filter(sampleID == "Modern") %>% 
  ggplot(aes(x=Length, y=Prop, colour = Genome)) + geom_line() + theme_bw() +
  xlab("Fragment length") + ylab("Proportion of reads") + ggtitle("Fragment length distribution for Modern")
#Plot for ModernL7, excluding M.oralis
normLengthDataQ30 %>% filter(sampleID == "Modern") %>% filter(Genome != "moralis") %>% 
  ggplot(aes(x=Length, y=Prop, colour = Genome)) + geom_line() + theme_bw() +
  xlab("Fragment length") + ylab("Proportion of reads") + 
  ggtitle("Fragment length distribution for Modern, excluding M.oralis")

###Plot distribution data as histogram, separating by genome

#Plot for Elsidron1
normLengthDataQ30 %>% filter(sampleID == "Elsidron1") %>% 
  ggplot(aes(x=Length, y=Prop, fill=Std)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_wrap((~Genome), ncol = 4) + theme(legend.position="none") + 
  scale_y_continuous("Proportion of total reads") +
  scale_x_continuous("Fragment length") +
  ggtitle("Fragment Lengths by Genome for Elsidron1")

#Plot for ModernL7
normLengthDataQ30 %>% filter(sampleID != "Elsidron1") %>% 
  ggplot(aes(x=Length, y=Prop, fill=Std)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_wrap((~Genome), ncol = 4) + theme(legend.position="none") + 
  scale_y_continuous("Proportion of total reads") +
  scale_x_continuous("Fragment length") +
  ggtitle("Fragment Lengths by Genome for Modern")

#Plot for Modern, excluding M.oralis
normLengthDataQ30 %>% filter(sampleID != "Elsidron1") %>% filter(Genome != "moralis") %>%
  ggplot(aes(x=Length, y=Prop, fill=Std)) + geom_bar(stat = "identity") + theme_bw() + 
  facet_wrap((~Genome), ncol = 4) + theme(legend.position="none") + 
  scale_y_continuous("Proportion of total reads") +
  scale_x_continuous("Fragment length") +
  ggtitle("Fragment Lengths by Genome for Modern, excluding M.oralis")


################# Plot overall length Data ##########################

#Box plot comparing overall fragment length data for each sample
lengthVectors_long %>% ggplot(aes(x=sampleID, y=lengths, fill=sampleID)) + stat_summary(fun.data = f, geom="boxplot") + 
  theme_bw() + theme(axis.title.x = element_blank()) + ylab("Average Fragment Length") + guides(fill=FALSE) +
  ggtitle("Average fragment lengths per sample")

#Distribution plot comparing fragment length data for each sample

#create object with summarised lengths for Elsidron1
elsLengthData <- lengthData %>% select(-FileName, -Genome, -Std) %>% filter(sampleID == "ELSIDRON1") %>% 
  select(-sampleID) %>% group_by(Length) %>% summarise_each(funs(sum)) %>% mutate(sampleID = "ELSIDRON1")
#Create object with summarised lengths for Modern sample
modLengthData <- lengthData %>% select(-FileName, -Genome, -Std) %>% filter(sampleID != "ELSIDRON1") %>% 
  select(-sampleID) %>% group_by(Length) %>% summarise_each(funs(sum)) %>% mutate(sampleID = "Modern")
#Bind objects into single data frame for plotting
lengthDist <- modLengthData %>% bind_rows(elsLengthData)
#Plot lengthDist
lengthDist %>% ggplot(aes(x=Length, y=Occurences, colour=sampleID)) + geom_line() + theme_bw() + 
  xlab("Read length") + ylab("Number of reads") + ggtitle("Distribution of fragment lengths for each sample")

################# Explain the reduction in mean fragment length for some genomes ##################

#Create a list of lengths for each genome, and determine the range in values (min & max)
modLgDataQ30Long <- lengthDataQ30Long %>% filter(sampleID == "Modern") %>% select(lengths, genome)
modLgRange <- split(modLgDataQ30Long, modLgDataQ30Long$genome) %>% lapply(function(x){range(x$lengths)}) %>% bind_rows()


######## plot lengths from fastq file #########

modFastqLengths <- read.csv(file="trimData/read_length.txt", sep="", header=FALSE)
names(modFastqLengths) <- c("occ", "length")
#Plot
modFastqLengths %>% ggplot(aes(x=length, y=occ)) + geom_line() + 
  theme_bw() + xlab("Fragment length") + ylab("Occurences")

#modFastqLengths <- modFastqLengths %>% mutate(freq = V1/58956478)
#modFastqLengths %>% sum(modFastqLengths$V1)
#[1] 58956478
#modFastqLenghts %>% ggplot(aes(x=V2, y=freq)) + geom_line() + theme_bw() + xlab("Fragment length") + ylab("Frequency")
#names(modFastqLenghts) <- c("occ", "length", "freq")
#modFastqFreq <- modFastqLenghts %>% select(length, freq)
#generate dataframe with aligned Modern data
#modBamLengths <- lengthDataQ30 %>% filter(sampleID == "Modern") %>% select(length, freq) %>% group_by(length) %>% 
#  summarise_each(funs(sum)) %>% mutate(freq = freq/1165054)
#modBamLengths <- as.data.frame(modBamLengths)
#modFastqBamLengthData <- left_join(modBamLengths, modFastqFreq, by = "length")
#names(modFastqBamLengthData) <- c("length", "bamFreq", "fastqFreq")

#Plot both data on top of each other.
#modFastqBamLengthData %>% ggplot(aes(x=length)) +
#  geom_line(aes(y=bamFreq, colour="bamFreq")) +
#  geom_line(aes(y=fastqFreq, colour="fastqFreq")) +
#  theme_bw() +
#  ylab("Frequency")

#Compare proportion of reads which aligned at different length ranges
mapLengths <- read.csv(file="alnData/mapLengths.txt", sep = "", header = FALSE) 
names(mapLengths) <- c("occ", "length")
rep(mapLengths$length, mapLengths$occ)

##Convert continuous data to categorical (i.e. summarise by length)
mapLengths$lengthRange <- cut(mapLengths$length, seq(21,191,10), right = FALSE, 
                              labels=c("21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", 
                                       "101-110", "111-120", "121-130", "131-140", "141-150", "151-160", "161-170", 
                                       "171-180", "181-190"))
#same for modFastqLengths
modFastqLengths$lengthRange <- cut(modFastqLengths$length, seq(21,191,10), right = FALSE, 
                                   labels=c("21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", 
                                            "101-110", "111-120", "121-130", "131-140", "141-150", "151-160", "161-170", 
                                            "171-180", "181-190"))

#Sum occurenced within each category
mapLengthsByRange <- mapLengths %>% select(-length) %>% group_by(lengthRange) %>% summarise_each(funs(sum))
fastqLengthsByRange <- modFastqLengths %>% select(-length) %>% group_by(lengthRange) %>% summarise_each(funs(sum))
#Bind both tibbles together
modLengthRanges <- left_join(mapLengthsByRange, fastqLengthsByRange, by = "lengthRange")
names(modLengthRanges) <- c("lengthRange", "bam", "fastq")
modLengthRanges <- modLengthRanges %>% mutate(propMapped = bam/fastq)
#Plot proportion mapped by length range
pdf("plots/prop.mapped.length.modern.pdf")
modLengthRanges %>% 
  ggplot(aes(x=lengthRange, y=propMapped)) + 
  geom_bar(stat = "identity", colour="#00BFC4", fill="#00BFC4") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  xlab("Length Range") +
  ylab("Proportion of total reads mapped") +
  ggtitle("Proportion of reads mapped, by length (Modern)")
dev.off()

################### Extract & plot high Qual length Data for 6 samples & 15 genomes ##########################

# create list of all .txt files in folder 
lgDistFiles <- list.files("mapData/highQualMapData/", pattern = "lgdistribution.txt", 
                          full.names = TRUE, recursive = TRUE)

#read-in data from each text file and bind into a data frame
lengthData <- lgDistFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>%
#    set_colnames(c("std", "length", "freq")) %>%
    mutate(fileName = x)
}) %>%
  bind_rows
lengthData <- lengthData %>% select(-X) 
names(lengthData) <- c("std", "length", "occ", "fileName")
#Remove unnecessary information from fileName
editLengthFileNames <- function(x){
  x$fileName <- gsub('mapData/highQualMapData//results_2NoAdapt_', '', x$fileName)
  x$fileName <- gsub('_split/lgdistribution.txt', '', x$fileName)
  x$fileName <- gsub('_150519_lACGTG_rATTGA', '', x$fileName)
  x$fileName <- gsub('L7_lTACTG_rCTCGA', '', x$fileName)
  x$fileName <- gsub('L7_lACTGT_rCTCGA', '', x$fileName)
  x$fileName <- gsub('L7_lAAGAG_rNONE', '', x$fileName)
  x$fileName <- gsub('L8_lGTACC_rCTCGA', '', x$fileName)
  x$fileName <- gsub('_L7', '', x$fileName)
}
#Apply function to lengthData
lengthData$fileName <- editLengthFileNames(lengthData)

#Split fileName into sampleID and genome
lengthData <- colsplit(lengthData$fileName, "_", names=c("sampleID", "genome")) %>% bind_cols(lengthData) %>% select(-fileName)

#specify labels for facet_wrap
labels <- c(CHIMP = "Chimp", ELSIDRON1 = "Elsidron1", ELSIDRON2 = "Elsidron2", Modern = "Modern", SPYNEW = "Spy II", SPYOLD = "Spy I")

#Plot lengths by genome (use facet_wrap for sampleID)
#pdf("plots/frag.length.genome.pdf")
lengthData %>% ggplot(aes(x=length, y=occ, colour=genome)) + 
  geom_line() + 
  theme_bw() + 
  facet_wrap(~sampleID, scale="free_y", ncol = 2, labeller = labeller(sampleID = labels)) +
  xlab("Fragment length") + 
  ylab("Frequency") + 
  ggtitle("Fragment lengths by genome (MAPQ > 30)")
#dev.off()

###Normalize and plot again

#Read-in text file and assign to object
highQsplitCount <- 
  read_delim(file="mapData/highQualMapData/high_qual_split_count.txt", delim = "\t", skip = 1, col_names = FALSE) %>% 
  set_colnames(c("fileName", "count"))

#Split fileName into sampleID and genome
#highQsplitCount <- colsplit(highQsplitCount$fileName, "_", names = c("Adapt", "sample")) %>% 
#  bind_cols(highQsplitCount) %>% select(-Adapt, -fileName)
#highQsplitCount <- colsplit(highQsplitCount$sample, "_l", names = c("sampleID", "extra")) %>%
#  bind_cols(highQsplitCount) %>% select(-sample)
#highQsplitCount <- colsplit(highQsplitCount$extra, "_r", names = c("misc1", "almostGenome")) %>%
#  bind_cols(highQsplitCount) %>% select(-extra, -misc1)
#highQsplitCount <- colsplit(highQsplitCount$almostGenome, "_", names = c("misc2", "genome")) %>%
#  bind_cols(highQsplitCount) %>% select(-almostGenome, -misc2)

#Use function editFileNames to edit the fileNames
source("editFileNames.R")
highQsplitCount$fileName <- editFileNames(highQsplitCount)
#Split edited fileName into sampleID and genome
highQsplitCount <- colsplit(highQsplitCount$fileName, "_", names = c("sampleID", "genome")) %>% bind_cols(highQsplitCount) %>% select(-fileName)

#Re-calculate into proportions and plot
highQsplitCount %>% select(-genome) %>% group_by(sampleID) %>% summarise_each(funs(sum)) -> totalQ30
names(totalQ30) <- c("sampleID", "totalCount")
highQsplitCount <- left_join(highQsplitCount, totalQ30, by = "sampleID")
#highQsplitCount %>% mutate(prop = count/totalCount) -> highQsplitCount

#First edit highQsplitCount so that sampleID identical to lengthData
#highQsplitCount$sampleID <- gsub('_150519', '', highQsplitCount$sampleID)
#highQsplitCount$sampleID <- gsub('L7', '', highQsplitCount$sampleID)
#highQsplitCount$sampleID <- gsub('L8', '', highQsplitCount$sampleID)
#highQsplitCount$sampleID <- gsub('_', '', highQsplitCount$sampleID)

#Split lengthData by sampleID and then by genome, summarise occ of each length and bind back together
lengthData %>% select(-std) %>% split(f = .$sampleID) %>% lapply(function(x){x %>% select(-sampleID) %>% split(f = .$genome) %>% 
      lapply(function(z){z %>% select(-genome) %>% group_by(length) %>% summarise_each(funs(sum))}) %>% bind_rows(.id = "genome")}) %>% 
  bind_rows(.id = "sampleID") -> lengthDataCollated

#bind total count column to lengthData
lengthDataCollated <- highQsplitCount %>% select(-totalCount) %>% left_join(lengthDataCollated)
#calculate proportion of lengths and plot
#pdf("plots/norm.fragment.length.genome.pdf")
lengthDataCollated %>% mutate(occ = occ/count) %>% ggplot(aes(x=length, y=occ, colour=genome)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~sampleID, ncol = 2, scales = "free", labeller = labeller(sampleID = labels)) + 
  xlab("Fragment length") + 
  ylab("Proportion of total reads") + 
  ggtitle("Normalized fragment lengths by genome (MAPQ > 30)")
#dev.off()

#Plot norm frag lengths for Elsidron1 alone
#pdf("plots/norm.frag.lengths.Elsidron1.pdf")
lengthDataCollated %>% filter(sampleID == "ELSIDRON1") %>% mutate(occ = occ/count) %>%
  ggplot(aes(x=length, y=occ, colour=genome)) +
  geom_line() +
  theme_bw() +
  xlab("Fragment length") +
  ylab("Proportion of reads") + 
  ggtitle("Normalized fragment lengths for Elsidron1 (MAPQ > 30)")
#dev.off()

#Plot norm frag lengths for Modern alone
#pdf("plots/norm.frag.lengths.Modern.pdf")
lengthDataCollated %>% filter(sampleID == "Modern") %>% mutate(occ = occ/count) %>%
  ggplot(aes(x=length, y=occ, colour=genome)) +
  geom_line() +
  theme_bw() +
  xlab("Fragment length") +
  ylab("Proportion of reads") + 
  ggtitle("Normalized fragment lengths for Modern (MAPQ > 30)")
#dev.off()

#Plot frag lengths each in genome Modern sample
#pdf("plots/compare.frag.lengths.genome.Modern.pdf")
lengthDataCollated %>% filter(sampleID == "Modern") %>% 
  ggplot(aes(x=length, y=occ, colour=genome)) + 
  geom_line() + theme_bw() + 
  facet_wrap(~genome, ncol=3, scales = "free") + 
  xlab("Fragment length") + 
  ylab("Number of reads") + 
  ggtitle("Compare fragment lengths of each species (Modern, MAPQ > 30)") + 
  guides(colour=FALSE) + 
  xlim(22,190)
#dev.off()

#Plot frag lengths each in genome Elsidron1 sample
#pdf("plots/compare.frag.lengths.genome.Elsidron1.pdf")
lengthDataCollated %>% filter(sampleID == "ELSIDRON1") %>% 
  ggplot(aes(x=length, y=occ, colour=genome)) + 
  geom_line() + theme_bw() + 
  facet_wrap(~genome, ncol=3, scales = "free") + 
  xlab("Fragment length") + 
  ylab("Number of reads") + 
  ggtitle("Compare fragment lengths by species (Elsidron1, MAPQ > 30)") + 
  guides(colour=FALSE) + 
  xlim(22,190)
#dev.off()


##Split data frame into separate samples
#lengthData %>% filter(sampleID == "CHIMP") %>% select(-sampleID) -> chimpLengthData
#lengthData %>% filter(sampleID == "ELSIDRON1") %>% select(-sampleID) -> elsidron1LengthData
#lengthData %>% filter(sampleID == "ELSIDRON2") %>% select(-sampleID) -> elsidron2LengthData
#lengthData %>% filter(sampleID == "Modern") %>% select(-sampleID) -> modernLengthData
#lengthData %>% filter(sampleID == "SPYOLD") %>% select(-sampleID) -> spy1LengthData

###For each new table need to remove strand, summarise lengths and finally calculate totals

#chimp
#chimpLengthData %>% select(-std) %>% 
#  split(f = .$genome) %>% 
#  lapply(function(x){x %>% 
#      select(length, occ) %>% 
#      group_by(length) %>% 
#      summarise_each(funs(sum))
#    }) %>% bind_rows(.id = "genome") -> chimpLengthData

#chimpLengthData %>% select(-length) %>% 
#  group_by(genome) %>% 
#  summarise_each(funs(sum)) %>% 
#  left_join(chimpLengthData, by = "genome") %>% 
#  mutate(occ.x = occ.y/occ.x) -> chimpLengthData

#names(chimpLengthData) <- c("genome", "prop", "length", "freq")

#elsidron1
#elsidron1LengthData %>% select(-std) %>% 
#  split(f = .$genome) %>% 
#  lapply(function(x){x %>% 
#      select(length, occ) %>% 
#      group_by(length) %>% 
#      summarise_each(funs(sum))
#  }) %>% bind_rows(.id = "genome") -> elsidron1LengthData

#elsidron1LengthData %>% select(-length) %>% 
#  group_by(genome) %>% 
#  summarise_each(funs(sum)) %>% 
#  left_join(elsidron1LengthData, by = "genome") %>% 
#  mutate(occ.x = occ.y/occ.x) -> elsidron1LengthData

#names(elsidron1LengthData) <- c("genome", "prop", "length", "freq")

#elsidron2
#elsidron2LengthData %>% select(-std) %>% 
#  split(f = .$genome) %>% 
#  lapply(function(x){x %>% 
#      select(length, occ) %>% 
#      group_by(length) %>% 
#      summarise_each(funs(sum))
#  }) %>% bind_rows(.id = "genome") -> elsidron2LengthData

#elsidron2LengthData %>% select(-length) %>% 
#  group_by(genome) %>% 
#  summarise_each(funs(sum)) 
#%>% 
#  left_join(elsidron2LengthData, by = "genome") %>% 
#  mutate(occ.x = occ.y/occ.x) -> elsidron2LengthData

#names(elsidron2LengthData) <- c("genome", "prop", "length", "freq")

#Modern
#modernLengthData %>% select(-std) %>% 
#  split(f = .$genome) %>% 
#  lapply(function(x){x %>% 
#      select(length, occ) %>% 
#      group_by(length) %>% 
#      summarise_each(funs(sum))
#  }) %>% bind_rows(.id = "genome") -> modernLengthData

#modernLengthData %>% select(-length) %>% 
#  group_by(genome) %>% 
#  summarise_each(funs(sum)) 
#%>% 
#  left_join(modernLengthData, by = "genome") %>% 
#  mutate(occ.x = occ.y/occ.x) -> modernLengthData

#names(modernLengthData) <- c("genome", "prop", "length", "freq")

#SpyI
#spy1LengthData %>% select(-std) %>% 
#  split(f = .$genome) %>% 
#  lapply(function(x){x %>% 
#      select(length, occ) %>% 
#      group_by(length) %>% 
#      summarise_each(funs(sum))
#  }) %>% bind_rows(.id = "genome") -> spy1LengthData

#spy1LengthData %>% select(-length) %>% 
#  group_by(genome) %>% 
#  summarise_each(funs(sum)) 
#%>% 
#  left_join(spy1LengthData, by = "genome") %>% 
#  mutate(occ.x = occ.y/occ.x) -> spy1LengthData

#names(spy1LengthData) <- c("genome", "prop", "length", "freq")

##Plot length proportions
#chimpLengthData %>% ggplot(aes(x=length, y=prop, colour=genome)) + geom_line() + theme_bw()
#elsidron1LengthData %>% ggplot(aes(x=length, y=prop, colour=genome)) + geom_line() + theme_bw()
#elsidron2LengthData %>% ggplot(aes(x=length, y=prop, colour=genome)) + geom_line() + theme_bw()
#modernLengthData %>% ggplot(aes(x=length, y=prop, colour=genome)) + geom_line() + theme_bw()
#spy1LengthData %>% ggplot(aes(x=length, y=prop, colour=genome)) + geom_line() + theme_bw()


######Box plots ############

##Chimp##
#use lapply to convert frequency data into vector of lengths for each genome
chimpLengthData %>% select(-prop) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$freq)}) -> chimpVector
#convert from list of vectors back into data frame
chimpVector <- names(chimpVector) %>% 
  lapply(function(x){data_frame(genome = x, lengths = chimpVector[[x]])}) %>% bind_rows()
#Add information about bacterial cell wall
cellWall <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                  "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                  "T.denticola", "T.forsythia"), 
                       cellWall = c("Gram +", "Gram +", "Gram -", "Gram +", "Gram -", "Gram -", "Mycobacterium",
                                    "Archeal", "Gram -", "Gram -", "Gram -", "Gram +", "Gram +", "Gram -", "Gram -"))
#Bind to chimpVector
chimpVector <- chimpVector %>% left_join(cellWall, by = "genome")

#Generate a function to apply to stat_summary so that can plot max and min values
quantiles <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#plot vectors as boxplot
#pdf("plots/chimp.length.boxplot.pdf")
chimpVector %>% 
  ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom="boxplot") +
  theme_bw() + ylab("Average fragment length") + xlab("") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  ggtitle("Chimp") +
  annotate("text", x = c(1:15), y=5, size=3,
          label = c("2,984", "239", "270", "78", "826", "18", "417", "51", "80", "2158", "166", "21", "27", "89", "3,884"))
#dev.off()

##Elsidron1##
#use lapply to convert frequency data into vector of lengths for each genome
elsidron1LengthData %>% select(-prop) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$freq)}) -> elsidron1Vector
#convert from list of vectors back into data frame
elsidron1Vector <- names(elsidron1Vector) %>% 
  lapply(function(x){data_frame(genome = x, lengths = elsidron1Vector[[x]])}) %>% bind_rows()
#Bind to chimpVector
elsidron1Vector <- elsidron1Vector %>% left_join(cellWall, by = "genome")
#Generate boxplot
#pdf("plots/elsidron1.length.boxplot.pdf")
elsidron1Vector %>%
  ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom="boxplot") +
  theme_bw() + ylab("Average fragment length") + xlab("") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  ggtitle("Elsidron1") +
  annotate("text", x = c(1:15), y=5, size=3,
           label = c("125,228", "1,367", "37,017", "8,344", "3,608", "279", "1,002", 
                     "87,895", "722", "7,113", "9,259", "3,005", "288", "11,343", "29,091"))
#dev.off()

##Elsidron2##
#use lapply to convert frequency data into vector of lengths for each genome
elsidron2LengthData %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> elsidron2Vector
#convert from list of vectors back into data frame
elsidron2Vector <- names(elsidron2Vector) %>% 
  lapply(function(x){data_frame(genome = x, lengths = elsidron2Vector[[x]])}) %>% bind_rows()
#Bind to chimpVector
elsidron2Vector <- elsidron2Vector %>% left_join(cellWall, by = "genome")
#Generate boxplot
#pdf("plots/elsidron2.length.boxplot.pdf")
elsidron2Vector %>%
  ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom="boxplot") +
  theme_bw() + ylab("Average fragment length") + xlab("") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  ggtitle("Elsidron2") +
  annotate("text", x = c(1:15), y=5, size=3,
           label = c("113,511", "125", "2,839", "86", "1,647", "102", "344", 
                     "2,440", "997", "648", "482", "945", "98", "1,864", "5,881"))
#dev.off()

##Modern##
modernLengthData %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> modernVector
#convert from list of vectors back into data frame
modernVector <- names(modernVector) %>% 
  lapply(function(x){data_frame(genome = x, lengths = modernVector[[x]])}) %>% bind_rows()
#Bind to chimpVector
modernVector <- modernVector %>% left_join(cellWall, by = "genome")
#Generate boxplot
#pdf("plots/modern.length.boxplot.pdf")
modernVector %>%
  ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom="boxplot") +
  theme_bw() + ylab("Average fragment length") + xlab("") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  ggtitle("Modern") +
  annotate("text", x = c(1:15), y=5, size=3,
           label = c("17,515", "534", "20,868", "155", "272,715", "1102", "967", 
                     "174", "32,278", "4,326", "44,039", "43,591", "421", "6,670", "86,641"))
#dev.off()

##SpyI##
spy1LengthData %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> spy1Vector
#convert from list of vectors back into data frame
spy1Vector <- names(spy1Vector) %>% 
  lapply(function(x){data_frame(genome = x, lengths = spy1Vector[[x]])}) %>% bind_rows()
#Bind to chimpVector
spy1Vector <- spy1Vector %>% left_join(cellWall, by = "genome")
#plot
pdf("plots/spy1.length.boxplot.pdf")
spy1Vector %>%
  ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom="boxplot") +
  theme_bw() + ylab("Average fragment length") + xlab("") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  theme(legend.position="bottom", legend.title = element_blank()) + 
  ggtitle("Spy I") +
  annotate("text", x = c(1:15), y=5, size=3,
           label = c("7,223", "75", "124", "46", "51", "74", "328", 
                     "1,487", "130", "51", "78", "875", "70", "57", "132"))
dev.off()

############################## Second attempt at code to create Box plots ###################################

#Split the lengthData data frame according to sampleID into separate data frames, each with name of the sample
list2env(split(lengthDataCollated, lengthDataCollated$sampleID), envir = .GlobalEnv)

##Attempts to create a function that will convert freq data to vector of lengths unsuccessful. 
#following works!
CHIMP %>% select(-sampleID, -count) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> CHIMP
CHIMP <- names(CHIMP) %>% 
  lapply(function(x){data_frame(genome = x, lengths = CHIMP[[x]])}) %>% bind_rows()

##Have to apply to each data frame individually

#Elsidron1
ELSIDRON1 %>% select(-sampleID, -count) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> ELSIDRON1
ELSIDRON1 <- names(ELSIDRON1) %>% 
  lapply(function(x){data_frame(genome = x, lengths = ELSIDRON1[[x]])}) %>% bind_rows()

#Elsidron2
ELSIDRON2 %>% select(-sampleID, -count) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> ELSIDRON2
ELSIDRON2 <- names(ELSIDRON2) %>% 
  lapply(function(x){data_frame(genome = x, lengths = ELSIDRON2[[x]])}) %>% bind_rows()

#SpyNew
na.omit(SPYNEW) -> SPYNEW
SPYNEW %>% select(-sampleID, -count) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> SPYNEW
SPYNEW <- names(SPYNEW) %>% 
  lapply(function(x){data_frame(genome = x, lengths = SPYNEW[[x]])}) %>% bind_rows()

#SpyOld
SPYOLD %>% select(-sampleID, -count) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> SPYOLD
SPYOLD <- names(SPYOLD) %>% 
  lapply(function(x){data_frame(genome = x, lengths = SPYOLD[[x]])}) %>% bind_rows()

#Modern
Modern %>% select(-sampleID, -count) %>% split(f = .$genome) %>% lapply(function(x){rep(x$length, x$occ)}) -> Modern
Modern <- names(Modern) %>% 
  lapply(function(x){data_frame(genome = x, length = Modern[[x]])}) %>% bind_rows()

#Generate vectors with information about cell wall type, GC content, oxidase tests, phylum
cellWall <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                  "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                  "T.denticola", "T.forsythia"), 
                       cellWall = c("GramPos", "GramPos", "GramNeg", "GramPos", "GramNeg", "GramNeg", "Mycobacterium",
                                    "Archeal", "GramNeg", "GramNeg", "GramNeg", "GramPos", "GramPos", "GramNeg", "GramNeg"))

GCcontent <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                   "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                   "T.denticola", "T.forsythia"), 
                        GC = c("65-70", "45-50", "45-50", "45-50", "<30", "35-40", "65-70", "<30", "50-55", "45-50", "40-45",
                               "40-45", "35-40", "35-40", "45-50"))

phylum <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                "T.denticola", "T.forsythia"), 
                     phylum = c("Actinobacteria", "Actinobacteria", "Proteobacteria", "Firmicutes", "Fusobacteria", "Proteobacteria", 
                                "Actinobacteria", "Euryarchaeota", "Proteobacteria", "Bacteroidetes", "Bacteroidetes", "Firmicutes", 
                                "Firmicutes", "Spirochaetes", "Bacteroidetes"))

catalase <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                 "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                 "T.denticola", "T.forsythia"), 
                      oxidase = c("V", "N", "N", "N", "N", "P", "P", "U", "P", "N", "N", "N", "N", "N", "V"))
                     
#Bind the above data to each data frame and then plot different variations
ELSIDRON1 %>% left_join(cellWall, by = "genome") %>% 
  left_join(phylum, by = "genome") %>% 
  left_join(GCcontent, by = "genome") %>% 
  left_join(catalase, by = "genome") -> ELSIDRON1

#Generate a function to apply to stat_summary so that can plot max and min values
quantiles <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#Add sampleSize data
Els1Counts <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                    "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                    "T.denticola", "T.forsythia"),
                         No.Reads = c("125228", "1367", "37017", "8344", "3608", "279", "1002", 
                                      "87895", "722", "7113", "9259", "3005", "288", "11343", "29091"))
ELSIDRON1 <- ELSIDRON1 %>% left_join(Els1Counts, by = "genome")
ELSIDRON1$No.Reads <- as.numeric(ELSIDRON1$No.Reads)

#Sort columns according to the variable wanting to plot             
ELSIDRON1 <- ELSIDRON1[order(ELSIDRON1$cellWall),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
ELSIDRON1$genome <- factor(ELSIDRON1$genome, levels = unique(ELSIDRON1$genome))
#Generate boxplot
ELSIDRON1 %>% ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1")

##### Repeat the plot, but remove genomes which have <1000 reads aligning
ELSIDRON1 %>% filter(No.Reads > 1000) %>%
  ggplot(aes(x=genome, y=lengths, fill=cellWall)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1 (> 1000 reads)")

#### New plot, compare according to phylum ###
ELSIDRON1 <- ELSIDRON1[order(ELSIDRON1$phylum),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
ELSIDRON1$genome <- factor(ELSIDRON1$genome, levels = unique(ELSIDRON1$genome))
ELSIDRON1 %>% ggplot(aes(x=genome, y=lengths, fill=phylum)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1")

#Repeat plot, but instead of using stat_summary, use boxplot
#ELSIDRON1 %>% ggplot(aes(x=genome, y=lengths, fill=phylum)) +
#  geom_boxplot(outlier.shape = NA) +
#  ylim(20,110) +
#  theme_bw() + 
#  ylab("Average fragment length") + xlab("") + 
#  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
#  ggtitle("Elsidron1")  

#Plot phylum, but remove genomes with low sample size
ELSIDRON1 %>% filter(No.Reads > 1000) %>%
  ggplot(aes(x=genome, y=lengths, fill=phylum)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1 (>1000 Reads)")
                     
#Plot according to oxidase activity
ELSIDRON1 <- ELSIDRON1[order(ELSIDRON1$oxidase),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
ELSIDRON1$genome <- factor(ELSIDRON1$genome, levels = unique(ELSIDRON1$genome))
#plot
ELSIDRON1 %>% ggplot(aes(x=genome, y=lengths, fill=oxidase)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1")
#Plot according to oxidase activity, Reads > 1000
ELSIDRON1 %>% filter(No.Reads > 1000) %>% 
  ggplot(aes(x=genome, y=lengths, fill=oxidase)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1 (>1000 Reads)")

#Plot according to GC content
ELSIDRON1 <- ELSIDRON1[order(ELSIDRON1$GC),]
#Convert genome to factor to prevent ggplot from re-ordering when plotting
ELSIDRON1$genome <- factor(ELSIDRON1$genome, levels = unique(ELSIDRON1$genome))
#plot
ELSIDRON1 %>% ggplot(aes(x=genome, y=lengths, fill=GC)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1")
#Replot with reads > 1000
ELSIDRON1 %>% filter(No.Reads > 1000) %>%
  ggplot(aes(x=genome, y=lengths, fill=GC)) + 
  stat_summary(fun.data = quantiles, geom = "boxplot") + 
  theme_bw() + 
  ylab("Average fragment length") + xlab("") + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  ggtitle("Elsidron1 (>1000 Reads)")

###Calculate key statistics about fragment lengths accroding to genome
Els1.lgStatsByGenome <- ELSIDRON1 %>% 
  select(genome, lengths) %>% 
  group_by(., genome) %>% 
  summarise(
    count = n(), 
    mean = mean(lengths, na.rm = TRUE), 
    sd = sd(lengths, na.rm = TRUE), 
    median = median(lengths, na.rm = TRUE), 
    IQR = IQR(lengths, na.rm = TRUE)
    )

##perform kruskal.wallis test
kruskal.test(lengths ~ genome, data = ELSIDRON1)
#Indicates there is a statistically significant difference in lengths.
#To determine which groups are statistically significant use pairwise.wilcox.test()
  #pairwise.wilcox.test(PlantGrowth$weight, PlantGrowth$group, p.adjust.method = "BH")
pairwise.wilcox.test(ELSIDRON1$lengths, ELSIDRON1$genome, p.adjust.method = "BH")

###Calculate key statistics about fragment lengths accroding to cellWally structure
Els1.lgStatsByCellWall <- ELSIDRON1 %>% 
  select(cellWall, lengths) %>% 
  group_by(., cellWall) %>% 
  summarise(
    count = n(), 
    mean = mean(lengths, na.rm = TRUE), 
    sd = sd(lengths, na.rm = TRUE), 
    median = median(lengths, na.rm = TRUE), 
    IQR = IQR(lengths, na.rm = TRUE)
  )




