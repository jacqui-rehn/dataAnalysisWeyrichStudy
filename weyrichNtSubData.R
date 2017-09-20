### Weyrich ntSub data

#load packages
library(dplyr)
library(readr)
library(magrittr)
library(tibble)
library(ggplot2)
library(reshape2)
library(stringr)

#Create a vector with desired column names
subDataColNames <- (c("Chr", "End", "Std", "Pos", "A", "C", "G", "T", "Total", 
                      "GtoA", "CtoT", "AtoG", "TtoC", "AtoC", "AtoT", "CtoG", 
                      "CtoA", "TtoG", "TtoA", "GtoC", "GtoT"))

# create list of all .txt files in folder 
ntSubFiles <- list.files("highQualAlnData", pattern = "misincorporation.txt",
                         full.names = TRUE, recursive = TRUE)

# read-in files and bind into a data frame that include the FileName
ntSubData <- ntSubFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>% select(1:21) %>%
    set_colnames(subDataColNames) %>% filter(Total > 0) %>% filter(Pos < 26) %>%
    mutate(FileName = x) %>% select(-Chr)
}) %>%
  bind_rows

#Edit FileName to include only the sampleID and Genome
ntSubData$FileName <- gsub('highQualAlnData/results_2NoAdapt_', '', ntSubData$FileName)
ntSubData$FileName <- gsub('_Q30_split/misincorporation.txt', '', ntSubData$FileName)
ntSubData$FileName <- gsub('ELSIDRON1', 'Elsidron1', ntSubData$FileName)

#Split FileName into separate columns - sampleID and Genome
ntSubData <- ntSubData %>% 
  mutate(Genome = str_extract(FileName, "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)"), 
         sampleID = str_replace(FileName, "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)", ""))

#Remove FileName column as no longer needed
ntSubData <- ntSubData %>% select(-FileName)
#Edit Genome names
ntSubData$Genome <- gsub('aoris', 'A.oris', ntSubData$Genome)
ntSubData$Genome <- gsub('cgracilis', 'C.gracilis', ntSubData$Genome)
ntSubData$Genome <- gsub('esaphenum', 'E.saphenum', ntSubData$Genome)
ntSubData$Genome <- gsub('fnucleatum', 'F.nucleatum', ntSubData$Genome)
ntSubData$Genome <- gsub('mneoaurum', 'M.neoaurum', ntSubData$Genome)
ntSubData$Genome <- gsub('pgingivalis', 'P.gingivalis', ntSubData$Genome)
ntSubData$Genome <- gsub('smutans', 'S.mutans', ntSubData$Genome)
ntSubData$Genome <- gsub('tdenticola', 'T.denticola', ntSubData$Genome)
ntSubData$Genome <- gsub('tforsythia', 'T.forsythia', ntSubData$Genome)
ntSubData$Genome <- gsub('moralis', 'M.oralis', ntSubData$Genome)

#Split data into separate samples and assign to own dataframes
Mod_ntSubData <- ntSubData %>% filter(sampleID == "Modern")
El1_ntSubData <- ntSubData %>% filter(sampleID == "Elsidron1")

########### Using a loop to plot data for all Elsidron1 samples  #############

# create graphing function
ntSub.graph <- function(df, na.rm = TRUE, ...){
  
  # Specify sampleID
  sampleID <- unique(df$sampleID)
  
  # create list of genomeID's in data to loop over 
  genomeID_list <- unique(df$Genome)
  
  # create for loop to split data based on sampleID 
  for (i in seq_along(genomeID_list)) {
      
      # create object to store 5p data
      SubFreq_5p <- subset(df, df$Genome==genomeID_list[i]) %>% filter(End == "5p") %>% 
        select(-End, -Std, -Genome, -sampleID) %>% 
        group_by(Pos) %>% summarise_each(funs(sum)) %>% 
        mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
               AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
               TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
        select(-A, -C, -G, -T, -Total) %>% 
        melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
      
      #create object to store 3p data
      SubFreq_3p <- subset(df, df$Genome==genomeID_list[i]) %>% filter(End == "3p") %>% 
        select(-End, -Std, -Genome, -sampleID) %>% 
        group_by(Pos) %>% summarise_each(funs(sum)) %>% 
        mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
               AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
               TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
        select(-A, -C, -G, -T, -Total) %>% 
        melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
      
      #plot to object, 5p data
      plot5pData <- SubFreq_5p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + geom_line() + 
        theme_bw() + scale_y_continuous(limits = c(0, 0.4), position = "left") + 
        ylab("Substitution Frequency") + xlab("Position from the 5' end") + theme(legend.position = "none")
      
      #plot to object, 3p data 
      plot3pData <- SubFreq_3p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + geom_line() + theme_bw() +
        theme(axis.title.y=element_blank()) + scale_x_reverse() + 
        scale_y_continuous(limits = c(0, 0.4), position = "right") +
        ylab("Substitution Frequency") + xlab("Position from the 3' end")
      
      #print plots
      grid.arrange(plot5pData, plot3pData, ncol=2, widths=c(0.85,1), 
                   top = (paste(genomeID_list[i], ' - nucleotide substitution data (', paste(sampleID, ')'))))
    
    #End loop
  }
}

# run graphing function on df
#pdf("plots/elsidron1.ntsub.genome.pdf")
ntSub.graph(El1_ntSubData) 
#dev.off()
#pdf("plots/modern.ntsub.genome.pdf")
ntSub.graph(Mod_ntSubData)
#dev.off()

################### Compare frequency of CtoT and GtoA mutations #############################

#Generate a list of 5pCtoT_freq.txt files
CtoT.Files <- list.files("highQualAlnData", pattern = "5pCtoT_freq.txt", full.names = TRUE, recursive = TRUE)
#Read in each of these files and bind together in a data frame
CtoTfreqData <- CtoT.Files %>% lapply(function(x){read_delim(x, delim = "\t", skip = 1, col_names = FALSE) %>% 
    set_colnames(c("Pos", "Freq")) %>% mutate(FileName = x)}) %>% bind_rows
#Edit FileName to include only sampleID and Genome
CtoTfreqData$FileName <- gsub('highQualAlnData/results_2NoAdapt_', '', CtoTfreqData$FileName)
CtoTfreqData$FileName <- gsub('_Q30_split/5pCtoT_freq.txt', '', CtoTfreqData$FileName)
#Split sampleID and Genome into separate vectors
CtoTfreqData <- CtoTfreqData %>% mutate(Genome = str_extract(FileName, 
                    "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|moralis|pgingivalis|smutans|tdenticola|tforsythia)"), 
                sampleID = str_replace(FileName, 
                    "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|moralis|pgingivalis|smutans|tdenticola|tforsythia)", "")) %>% 
  select(-FileName)
#Add information on cell.wall and bind to data frame
cell.wall <- data_frame(Genome = c("aoris", "cgracilis", "esaphenum", "fnucleatum", 
                                   "mneoaurum","moralis", "pgingivalis", "smutans", 
                                   "tdenticola", "tforsythia"), 
                        cell.wall = c("gram+", "gram-", "gram+", "gram-", "Mycobacterium",
                                      "Archeal", "gram-", "gram+", "gram-", "gram-"))
CtoTfreqData <- CtoTfreqData %>% left_join(cell.wall, by = "Genome")

#Generate a dot plot representing the CtoT frequency at position 1
pdf("plots/CtoTfreq.mod.els1.pdf")
CtoTfreqData %>% filter(Pos == "1") %>% ggplot(aes(x=Genome, y=Freq, shape=sampleID, colour=cell.wall)) + 
  geom_point(size=4) + scale_shape_manual(values=c(15,13), labels=c("Elsidron1", "Modern")) + 
  theme_bw() + ylab("Misincorporation Frequency") + 
  ggtitle("Cytosine to Thymine misincorporation frequency\nat Position 1 of the 5' end") + 
  scale_x_discrete(labels=c("A.oris", "C.gracilis", "E.saphenum", "F.nucleatum", "M.neoaurum", "M.oralis", 
                            "P.gingivalis", "S.mutans", "T.denticola", "T.forsythia")) + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + labs(shape="Sample ID", colour="Cell Wall") + 
  scale_colour_discrete(labels=c("Archeal", "Gram Positive", "Gram Negative", "Mycobacterium"))
dev.off()

#Plot average of last 5 positions
xlabs <- c("A.oris", "C.gracilis", "E.saphenum", "F.nucleatum", "M.neoaurum", "M.oralis", 
           "P.gingivalis", "S.mutans", "T.denticola", "T.forsythia")
#Write a function for generating boxplot
singleNtFreqBoxPlot <- function(x){x %>% ggplot(aes(x=Genome, y=Freq, fill=cell.wall)) +
    geom_boxplot() + facet_wrap(~sampleID, ncol=1) + theme_bw() + ylab("Misincorporation Frequency") + ylim(0,0.55) +
    scale_x_discrete(labels=(xlabs)) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) +
    scale_fill_discrete(labels=c("Archeal", "Gram Positive", "Gram Negative", "Mycobacterium")) + labs(fill="Cell Wall")}

CtoTfreqData %>% filter(Pos < 4) %>% ggplot(aes(x=Genome, y=Freq, fill=cell.wall)) + 
  geom_boxplot() + facet_wrap(~sampleID, ncol=1) + theme_bw() + 
  ylab("Misincorporation Frequency") + ggtitle("C->T misincorporation frequency, averaged of 3 bases at 5' end") + ylim(0,0.55) + 
  scale_x_discrete(labels=c(xlabs)) + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=1)) +
  scale_fill_discrete(labels=c("Archeal", "Gram Positive", "Gram Negative", "Mycobacterium")) + labs(fill="Cell Wall")

###Generate sample plots but with G->A mutation at 3p end

#Generate a list of 5pCtoT_freq.txt files
GtoA.files <- list.files("highQualAlnData", pattern = "3pGtoA_freq.txt", full.names = TRUE, recursive = TRUE)
#Read in data from each text file and bind into a single data frame
##Most values in Modern_moralis_Q30 were 0 which affected ability to read-in data (values not 0 were converted to NA)
##Solved by specifying col_types as number
GtoAfreqData <- GtoA.files %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 1, col_names = FALSE, col_types = cols("i", "n")) %>% 
    set_colnames(c("Pos", "Freq")) %>% mutate(FileName = x)
  }) %>% bind_rows
#Edit FileName to include only sampleID and Genome
GtoAfreqData$FileName <- gsub('highQualAlnData/results_2NoAdapt_', '', GtoAfreqData$FileName)
GtoAfreqData$FileName <- gsub('_Q30_split/3pGtoA_freq.txt', '', GtoAfreqData$FileName)

#Write function to split FileName into separate columns - sampleID and Genome
splitFileName <- function(x){mutate(x, Genome = str_extract(FileName, 
          "(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)"), 
          sampleID = str_replace(FileName, "_(aoris|cgracilis|esaphenum|fnucleatum|mneoaurum|pgingivalis|smutans|tdenticola|tforsythia|moralis)", ""))}
#Apply to GtoAfreqData
GtoAfreqData <- splitFileName(GtoAfreqData)
#Add cell wall information with left_join
GtoAfreqData <- GtoAfreqData %>% left_join(cell.wall, by = "Genome")

#Plot misincorporation frequency at Position 1
pdf("plots/GtoAfreq.mod.els1.pdf")
GtoAfreqData %>% filter(Pos == "1") %>% ggplot(aes(x=Genome, y=Freq, shape=sampleID, colour=cell.wall)) + 
  geom_point(size=4) + theme_bw() + scale_shape_manual(values=c(15,13), labels=c("Elsidron1", "Modern")) + 
  ylab("Misincorporation Frequency") + ggtitle("Guanine->Adenine misincorporation frequency\nat position 1 of the 3' end") + 
  scale_x_discrete(labels=(xlabs)) + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 1)) + 
  labs(shape="Sample ID", colour="Cell Wall") + 
  scale_color_discrete(labels=c("Archeal", "Gram Positive", "Gram Negative", "Mycobacterium"))
dev.off()

#Plot average misincorporation frequency of the last 3 positions
GtoAfreqData %>% filter(Pos < 4) %>% singleNtFreqBoxPlot() + ggtitle("G->A misincorporation frequency, averaged of 3 bases at 3' end")

###################Compare overall C->T and G->A frequency for each sample#################

###Summarise Frequency of G->A mutations over last 10 bases at 3' end
#Summarise data for Elsidron1 and assign to object
GtoA.Elsidron <- GtoAfreqData %>% filter(Pos < 11) %>% select(-FileName, -Genome, -cell.wall) %>% 
  filter(sampleID == "ELSIDRON1") %>% select(-sampleID) %>% group_by(Pos) %>% summarise_each(funs(sum)) %>% 
  mutate(sampleID = "Elsidron1")
#Summarise data for Modern and assign to object
GtoA.Modern <- GtoAfreqData %>% filter(Pos < 11) %>% select(-FileName, -Genome, -cell.wall) %>% 
  filter(sampleID == "Modern") %>% select(-sampleID) %>% group_by(Pos) %>% summarise_each(funs(sum)) %>%
  mutate(sampleID = "Modern")
#Bind both objects into a single data frame
overallGtoA <- GtoA.Elsidron %>% bind_rows(GtoA.Modern)

###Summarise Frequency of C->T mutations over first 10 bases at 5' end
#Summarise data for Elsidron1
CtoT.Elsidron <- CtoTfreqData %>% filter(Pos < 11) %>% select(-FileName, -Genome, -cell.wall) %>% 
  filter(sampleID == "ELSIDRON1") %>% select(-sampleID) %>% group_by(Pos) %>% summarise_each(funs(sum)) %>% 
  mutate(sampleID = "Elsidron1")
#Summarise data for Modern
CtoT.Modern <- CtoTfreqData %>% filter(Pos < 11) %>% select(-FileName, -Genome, -cell.wall) %>% 
  filter(sampleID == "Modern") %>% select(-sampleID) %>% group_by(Pos) %>% summarise_each(funs(sum)) %>% 
  mutate(sampleID = "Modern")
#Bind into a single data frame
overallCtoT <- CtoT.Elsidron %>% bind_rows(CtoT.Modern)

#Plot overall data for both mutations and arrange side-by-side

CtoTPlot <- overallCtoT %>% ggplot(aes(x=Pos, y=Freq, colour=sampleID)) + geom_point() + theme_classic() + 
  stat_smooth(se = FALSE) + xlab("Position") + ylab("Misincorporation Frequency") + 
  ggtitle("C->T misincorporation at 5' end") + ylim(0, 4) + guides(colour=FALSE)

GtoAPlot <- overallGtoA %>% ggplot(aes(x=Pos, y=Freq, colour=sampleID)) + geom_point() + theme_classic() + 
  stat_smooth(se = FALSE) + scale_x_reverse() + scale_y_continuous(limits = c(0, 4), position = "right") + 
  ggtitle("G->A misincorporation at 3' end") + xlab("Position") + ylab("Misincorporation Frequency") + 
  theme(legend.position = c(0.2,0.8)) + labs(colour="Sample ID")

library(gridExtra)
grid.arrange(CtoTPlot, GtoAPlot, ncol=2)


##################### New plots for expanded sample/species run ############################

#Create a vector with desired column names
subDataColNames <- (c("Chr", "End", "Std", "Pos", "A", "C", "G", "T", "Total", 
                      "GtoA", "CtoT", "AtoG", "TtoC", "AtoC", "AtoT", "CtoG", 
                      "CtoA", "TtoG", "TtoA", "GtoC", "GtoT"))

# create list of all .txt files in folder 
ntSubFiles <- list.files("mapData/highQualMapData/", pattern = "misincorporation.txt",
                         full.names = TRUE, recursive = TRUE)
# exclude SpyNEW as too few reads to analyse damage patterns
ntSubFiles <- ntSubFiles[!grepl("SPYNEWL8", ntSubFiles)]

# read-in files and bind into a data frame that includes the fileName
ntSubData <- ntSubFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE) %>% select(1:21) %>%
    set_colnames(subDataColNames) %>% filter(Total > 0) %>% filter(Pos < 26) %>%
    mutate(fileName = x) %>% select(-Chr)
}) %>%
  bind_rows

#Edit fileName to include sampleID and Genome
#colsplit(ntSubData, "_", names = c("directory", "adapters", "sampleID", "other"))
ntSubData$fileName <- gsub("mapData/highQualMapData//results_2NoAdapt_", '', ntSubData$fileName)
ntSubData$fileName <- gsub('_split/misincorporation.txt', '', ntSubData$fileName)
ntSubData$fileName <- gsub('_150519_lACGTG_rATTGA', '', ntSubData$fileName)
ntSubData$fileName <- gsub('L7_lACTGT_rCTCGA', '', ntSubData$fileName)
ntSubData$fileName <- gsub('L7_lTACTG_rCTCGA', '', ntSubData$fileName)
ntSubData$fileName <- gsub('L7_lAAGAG_rNONE', '', ntSubData$fileName)
ntSubData$fileName <- gsub('_L7L8_lGTACC_rCTCGA', '', ntSubData$fileName)

#Split fileName into sampleID and genome
colsplit(ntSubData$fileName, "_", names = c("sampleID", "genome")) %>% bind_cols(ntSubData) %>% select(-fileName) -> ntSubData

#Split data into separate samples and assign to own dataframes
Mod_ntSubData <- ntSubData %>% filter(sampleID == "Modern")
El1_ntSubData <- ntSubData %>% filter(sampleID == "ELSIDRON1")
El2_ntSubData <- ntSubData %>% filter(sampleID == "ELSIDRON2")
Chimp_ntSubData <- ntSubData %>% filter(sampleID == "CHIMP")
SpyOld_ntSubData <- ntSubData %>% filter(sampleID == "SPYOLD")

# create graphing function
ntSub.graph <- function(df, na.rm = TRUE, ...){
  
  # Specify sampleID
  sampleID <- unique(df$sampleID)
  
  # create list of genomeID's in data to loop over 
  genomeID_list <- unique(df$genome)
  
  # create for loop to split data based on sampleID 
  for (i in seq_along(genomeID_list)) {
    
    # create object to store 5p data
    SubFreq_5p <- subset(df, df$genome==genomeID_list[i]) %>% filter(End == "5p") %>% 
      select(-End, -Std, -genome, -sampleID) %>% 
      group_by(Pos) %>% summarise_each(funs(sum)) %>% 
      mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
             AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
             TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
      select(-A, -C, -G, -T, -Total) %>% 
      melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
    
    #create object to store 3p data
    SubFreq_3p <- subset(df, df$genome==genomeID_list[i]) %>% filter(End == "3p") %>% 
      select(-End, -Std, -genome, -sampleID) %>% 
      group_by(Pos) %>% summarise_each(funs(sum)) %>% 
      mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
             AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
             TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
      select(-A, -C, -G, -T, -Total) %>% 
      melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
    
    #plot to object, 5p data
    plot5pData <- SubFreq_5p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + geom_line() + 
      theme_bw() + scale_y_continuous(limits = c(0, 0.4), position = "left") + 
      ylab("Substitution Frequency") + xlab("Position from the 5' end") + theme(legend.position = "none")
    
    #plot to object, 3p data 
    plot3pData <- SubFreq_3p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + geom_line() + theme_bw() +
      theme(axis.title.y=element_blank()) + scale_x_reverse() + 
      scale_y_continuous(limits = c(0, 0.4), position = "right") +
      ylab("Substitution Frequency") + xlab("Position from the 3' end")
    
    #print plots
    library(gridExtra)
    grid.arrange(plot5pData, plot3pData, ncol=2, widths=c(0.85,1), 
                 top = (paste(genomeID_list[i], ' - nucleotide substitution data (', paste(sampleID, ')'))))
    
    #End loop
  }
}


# run graphing function on df
#pdf("plots/elsidron1.ntsub.genome.expanded.pdf")
ntSub.graph(El1_ntSubData) 
#dev.off()

#pdf("plots/modern.ntsub.genome.expanded.pdf")
ntSub.graph(Mod_ntSubData)
#dev.off()

#pdf("plots/elsidron2.ntsub.genome.pdf")
ntSub.graph(El2_ntSubData)
#dev.off()

#pdf("plots/SpyOld.ntsub.genome.pdf")
ntSub.graph(SpyOld_ntSubData)
#dev.off()

#pdf("plots/chimp.ntsub.genome.pdf")
ntSub.graph(Chimp_ntSubData)
#dev.off()

############### mutation frequency comparison between genomes of each sample ####################

#Generate a list of 5pCtoT_freq.txt files
CtoT.Files <- list.files("mapData/highQualMapData", pattern = "5pCtoT_freq.txt", full.names = TRUE, recursive = TRUE)
# exclude SpyNEW as too few reads to analyse damage patterns
CtoT.Files <- CtoT.Files[!grepl("SPYNEWL8", CtoT.Files)]
#Read in each of these files and bind together in a data frame (MUST specify col_types!)
CtoTfreqData <- CtoT.Files %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 1, col_names = FALSE, col_types = cols("i", "n")) %>% 
    set_colnames(c("pos", "freq")) %>% 
    mutate(fileName = x)}) %>% 
  bind_rows()

#Edit fileName to include only sampleID and genome
CtoTfreqData$fileName <- gsub('mapData/highQualMapData/results_2NoAdapt_', '', CtoTfreqData$fileName)
CtoTfreqData$fileName <- gsub('_split/5pCtoT_freq.txt', '', CtoTfreqData$fileName)
CtoTfreqData$fileName <- gsub('_150519_lACGTG_rATTGA', '', CtoTfreqData$fileName)
CtoTfreqData$fileName <- gsub('L7_lACTGT_rCTCGA', '', CtoTfreqData$fileName)
CtoTfreqData$fileName <- gsub('L7_lTACTG_rCTCGA', '', CtoTfreqData$fileName)
CtoTfreqData$fileName <- gsub('L7_lAAGAG_rNONE', '', CtoTfreqData$fileName)
CtoTfreqData$fileName <- gsub('_L7L8_lGTACC_rCTCGA', '', CtoTfreqData$fileName)
#Split sampleID and Genome into separate vector
colsplit(CtoTfreqData$fileName, "_", names = c("sampleID", "genome")) %>% bind_cols(CtoTfreqData) %>% select(-fileName) -> CtoTfreqData

#Add information on cell.wall and bind to data frame
#Create data-frame with information about bacterial cell wall
cellWall <- data_frame(genome = c("A.oris", "A.parvulum", "C.gracilis", "E.saphenum", "F.nucleatum", "H.influenza",
                                  "M.neoaurum","M.oralis", "N.meningitidis", "P.gingivalis", "P.intermedia", "S.mitis", "S.mutans", 
                                  "T.denticola", "T.forsythia"), 
                       cellWall = c("Gram +", "Gram +", "Gram -", "Gram +", "Gram -", "Gram -", "Mycobacterium",
                                    "Archeal", "Gram -", "Gram -", "Gram -", "Gram +", "Gram +", "Gram -", "Gram -"))
#Bind this information to substitution frequency data
CtoTfreqData <- CtoTfreqData %>% left_join(cellWall, by = "genome")
#Sort data based on sampleID [,1] then cellWall [,5]
CtoTfreqData <- CtoTfreqData[order( CtoTfreqData[,1], CtoTfreqData[,5] ), ]

#Convert to factor to prevent ggplot from re-ordering when plotting (Will get a warning message)
#CtoTfreqData$cellWall <- factor(CtoTfreqData$cellWall, levels = CtoTfreqData$cellWall)
#CtoTfreqData$cellWall <- factor(CtoTfreqData$cellWall, levels=unique(CtoTfreqData$cellWall))

#Convert genome to factor to prevent ggplot from re-ordering when plotting
CtoTfreqData$genome <- factor(CtoTfreqData$genome, levels = unique(CtoTfreqData$genome))

#Subset Elsidron1 data and plot
#pdf("plots/CtoT.Elsidron1.pdf")
CtoTfreqData %>% 
  filter(sampleID == "ELSIDRON1") %>% 
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Cytosine to Thymine misincorporation frequency\nat Position 1 of the 5' end") + 
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 0.255, ymax = 0.27, alpha = .2) + 
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = 0.16, ymax = 0.18, alpha = .2)
#dev.off()

#Plot all data on same axes
#plotCtoT <- 
#pdf("plots/CtoTfreq.overall.pdf")
CtoTfreqData %>% 
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ylim(0, 0.52) + 
#  theme(legend.position = "none") +
  ggtitle("Cytosine to Thymine misincorporation frequency\nat Position 1 of the 5' end")
#dev.off()

#Generate a list of 3pGtoA_freq.txt files
GtoA.Files <- list.files("mapData/highQualMapData", pattern = "3pGtoA_freq.txt", full.names = TRUE, recursive = TRUE)
# exclude SpyNEW as too few reads to analyse damage patterns
GtoA.Files <- GtoA.Files[!grepl("SPYNEWL8", GtoA.Files)]
#Read in each of these files and bind together in a data frame (MUST specify col_types!)
GtoAfreqData <- GtoA.Files %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 1, col_names = FALSE, col_types = cols("i", "n")) %>% 
    set_colnames(c("pos", "freq")) %>% 
    mutate(fileName = x)}) %>% 
  bind_rows()

#Edit fileName to include only sampleID and genome
GtoAfreqData$fileName <- gsub('mapData/highQualMapData/results_2NoAdapt_', '', GtoAfreqData$fileName)
GtoAfreqData$fileName <- gsub('_split/3pGtoA_freq.txt', '', GtoAfreqData$fileName)
GtoAfreqData$fileName <- gsub('_150519_lACGTG_rATTGA', '', GtoAfreqData$fileName)
GtoAfreqData$fileName <- gsub('L7_lACTGT_rCTCGA', '', GtoAfreqData$fileName)
GtoAfreqData$fileName <- gsub('L7_lTACTG_rCTCGA', '', GtoAfreqData$fileName)
GtoAfreqData$fileName <- gsub('L7_lAAGAG_rNONE', '', GtoAfreqData$fileName)
GtoAfreqData$fileName <- gsub('_L7L8_lGTACC_rCTCGA', '', GtoAfreqData$fileName)
#Split sampleID and Genome into separate vector
colsplit(GtoAfreqData$fileName, "_", names = c("sampleID", "genome")) %>% bind_cols(GtoAfreqData) %>% select(-fileName) -> GtoAfreqData

#Add information on cell.wall and bind to data frame
GtoAfreqData <- GtoAfreqData %>% left_join(cellWall, by = "genome")
#Sort data based on sampleID [,1] then cellWall [,5]
GtoAfreqData <- GtoAfreqData[order( GtoAfreqData[,1], GtoAfreqData[,5] ), ]

#Convert genome to factor to prevent ggplot from re-ordering when plotting
GtoAfreqData$genome <- factor(GtoAfreqData$genome, levels = unique(GtoAfreqData$genome))

#Plot all data on same axes
#plotGtoA <- 
#pdf("plots/GtoA.freq.overall.pdf")
GtoAfreqData %>% 
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  ylab("Misincorporation frequency") + 
  ylim(0,0.52) +
  ggtitle("Guanine to Adenine misincorporation frequency\nat Position 1 of the 3' end") 
#dev.off()

##plot 5' C->T next to 3' G->A
#grid.arrange(plotCtoT, plotGtoA, ncol=2, widths=c(0.8,1))


############### Add info on sample size ##################

#Read-in count data
highQsplitCount <- read_delim("mapData/highQualMapData/high_qual_split_count.txt", 
                              delim = "\t", skip = 1, col_names = FALSE) %>% 
  set_colnames(c("fileName", "count"))

###Edit fileName and split into sampleID and genome
#Edit fileName
highQsplitCount$fileName <- gsub('2NoAdapt_', '', highQsplitCount$fileName)
highQsplitCount$fileName <- gsub('_150519_lACGTG_rATTGA', '', highQsplitCount$fileName)
highQsplitCount$fileName <- gsub('L7_lACTGT_rCTCGA', '', highQsplitCount$fileName)
highQsplitCount$fileName <- gsub('L7_lTACTG_rCTCGA', '', highQsplitCount$fileName)
highQsplitCount$fileName <- gsub('L7_lAAGAG_rNONE', '', highQsplitCount$fileName)
highQsplitCount$fileName <- gsub('_L7L8_lGTACC_rCTCGA', '', highQsplitCount$fileName)
highQsplitCount$fileName <- gsub('L8_lGTACC_rCTCGA', '', highQsplitCount$fileName)
#Split sampleID and Genome into separate vector
colsplit(highQsplitCount$fileName, "_", names = c("sampleID", "genome")) %>% bind_cols(highQsplitCount) %>% select(-fileName) -> highQsplitCount
#exclude SPYNEW
highQsplitCount %>% filter(sampleID != "SPYNEW") -> highQsplitCount

#Bind sample counts to both the CtoT and GtoA ntSubData **WARNING - left_join will coerce genome from factor back to character!!
CtoTfreqData <- left_join(CtoTfreqData, highQsplitCount)
GtoAfreqData <- left_join(GtoAfreqData, highQsplitCount)

#Convert genome back to factor to prevent ggplot from re-ordering when plotting
GtoAfreqData$genome <- factor(GtoAfreqData$genome, levels = unique(GtoAfreqData$genome))
CtoTfreqData$genome <- factor(CtoTfreqData$genome, levels = unique(CtoTfreqData$genome))

#Filter out results where the sample size is < 100bp and plot:
pdf("plots/CtoTfreq.100more.reads.pdf")
CtoTfreqData %>% 
  filter(count > 100) %>% 
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Cytosine to thymine misincorporation frequency\n(> 100 reads per genome)")
dev.off()

pdf("plots/GtoAfreq.100more.reads.pdf")
GtoAfreqData %>% 
  filter(count > 100) %>% 
  filter(pos == "1") %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Guanine to adenine misincorporation frequency\n(>100 reads per genome)")
dev.off()

#Filter out results where the sample size is <500bp and plot:
pdf("plots/CtoTfreq.500more.reads.pdf")
CtoTfreqData %>% 
  filter(pos == "1") %>% 
  filter(count > 500) %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Cytosine to thymine misincorporation frequency\n(>500 reads per genome)")
dev.off()

pdf("plots/GtoAfreq.500more.reads.pdf")
GtoAfreqData %>% 
  filter(pos == "1") %>% 
  filter(count > 500) %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Guanine to adenine misincorporation frequency\n(>500 reads per genome)")
dev.off()

#Filter out results where the sample size is <1000 bp and plot:
CtoTfreqData %>% 
  filter(pos == "1") %>% 
  filter(count > 1000) %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Cytosine to thymine misincorporation frequency\n (>1000 reads per genome)")

GtoAfreqData %>% 
  filter(pos == "1") %>% 
  filter(count > 1000) %>% 
  ggplot(aes(x=genome, y=freq, shape=sampleID, colour=cellWall)) + 
  geom_point(size=4) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  ylab("Misincorporation frequency") + 
  ggtitle("Guanine to adenine misincorporation frequency\n (>1000 reads per genome)")
