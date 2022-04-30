# import seq data
seqData <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_consensus (1).xlsx")

# filter for insect OTUs
seqInsecta <- which(seqData$class=="Insecta")
insectData <- seq.data[c(seqInsecta),]

# remove unwanted columns
insectData <- insectData[,-c(56,185,248,249,272,274:275,277:291)]

library(writexl)
writexl::write_xlsx(insectData, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\insect.filtered.data.xlsx")

# wPOO function
#####
wPOO_function <- function(df) {                                                                              
  for (row in 1:nrow(df)) {                                                                
    for (col in 1:ncol(df)) {                                                              
      if (df[row, col] > 0) {                                                            
        df[row, col] <- 1
      }                                                                                                    
    }                                                                                                      
  }
  
  insectSampleDiet <- c()                                                                                
  for (seq_col in 1:ncol(df)) {                                                               
    J <- 0
    for (seq_row in 1:nrow(df)) {                                                             
      ifelse(df[seq_row, seq_col] > 0, J <- J + 1, J <- J)                                    
    }                                                                                                      
    insectSampleDiet <- append(insectSampleDiet, J)                                                    
  }
  colnames(df) <- c(insectSampleDiet)                                                       
  
  df <- as.data.frame(df)                                                        
  reads <- c()
  for (insect_col in 1:ncol(df)) {                                                            
    reads <- append(reads, df[,insect_col])                                                   
  }                                                                                                        
  
  insectDF <- data.frame(Row = rep(1:length(df[,1]), times = length(df[1,])),                                                   
                          Reads = reads,                                                                   
                          TotalReads = rep(insectSampleDiet, each = length(df[,1])))                               
  
  insectDF <- cbind(insectDF, insectDF$Reads/insectDF$TotalReads)                                     
  
  insectValueDF <- df                                                                       
  for (row in 1:nrow(insectValueDF)) {
    new_insectDF <- insectDF[insectDF$Row==row,]                                                        
    insectValueDF[row,] <- new_insectDF$`insectDF$Reads/insectDF$TotalReads`                         
  }                                                                                                        
  
  insect_wPOO <- c()                                                                                       
  for (value_row in 1:nrow(insectValueDF)) {                                                             
    insect_wPOO <- append(insect_wPOO, 100*(1/length(insectSampleDiet))*sum(insectValueDF[value_row,]))
  }  
  return(insect_wPOO)
}
#####

# import habitat data
# select columns of interest and remove NAs
#####
farms <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\FecesSampleDatabase6Nov2019_MISEQ1&MISEQ2.xlsx", sheet = "Samples")                        
farms <- data.frame(farms$Lab.nbr...26, farms$animal, farms$Species, farms$Location,
                    farms$Site, as.character.POSIXt(farms$Date))                     
farms <- farms[-145,]                                                                

farms_na <- which(is.na(farms$farms.Lab.nbr...26), arr.ind=TRUE)                     
farms <- farms[-c(farms_na),]

sort_farms <- farms[c(4,5)]                                                          
farm_names <- unique(sort_farms$farms.Location)                                      
farm_names <- farm_names[!is.na(farm_names)]
#####

# use hash match up missing farm names
#####
library(hash)
nameMap <- hash()                                

for (r in 1:nrow(sort_farms)) {                  
  locationName <- sort_farms[r, "farms.Location"]
  siteName <- sort_farms[r, "farms.Site"]        
  nameMap[[siteName]] <- locationName            
}
#####

# find lab numbers common to both datasets to identify which samples were actually sequenced
#####
farmsLabNbrs <- as.character(farms$farms.Lab.nbr...26)
zbj_names <- names(insectData[,2:268])           

LabNbr <- c()                                     
for (n in zbj_names) {                            
  for (l in farmsLabNbrs) {                           
    if (l==n)
      LabNbr <- append(LabNbr, as.numeric(n))     
  }                                               
}
#####

# compile information of interest from both datasets
#####
library(stringr)
zbj_df <- data.frame()                                                                                                                  

for (i in LabNbr) {                                                                                                                     
  seqSamples <- farms[which(farms$farms.Lab.nbr...26==i),]                                                                              
  zbj_df <- rbind(zbj_df, seqSamples) 
}    

for (d in 1:nrow(zbj_df)) {                                                                                                             
  if (is.na(zbj_df[d, "farms.Location"])) {                                                                                             
    zbj_df[d, "farms.Location"] <- nameMap[[zbj_df[d, "farms.Site"]]]                                                                   
  }
  zbj_df[d, "as.character.POSIXt.farms.Date."] <- ifelse(str_detect(zbj_df[d, "as.character.POSIXt.farms.Date."], "2017"), "wet", "dry")
}

zbj_final <- zbj_df[order(zbj_df[,4]),]                                                                                                 
colnames(zbj_final) <- c("Lab.nbr", "Animal", "Species", "Location", "Site", "Season")
#####

write_xlsx(zbj_final, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_final.xlsx")

insectData <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\insect.filtered.data.xlsx")
zbj_final <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_final.xlsx")

# create hash for matching lab numbers to species name
#####
predator_names <- hash()                             
for (row in 1:nrow(zbj_final)) {                     
  predNbr <- as.character(zbj_final[row, "Lab.nbr"])
  predSp <- as.character(zbj_final[row,"Species"])  
  predator_names[[predNbr]] <- predSp              
}
#####

# use hash to create a matrix of sequence data for the focal species
#####
zbjInsectData <- insectData[,c(2:268)]
zbjInsectData <- zbjInsectData[,-c(which(names(zbjInsectData) %in% c("146","148")))]

# remove unwanted columns
new_zbjInsectData <- zbjInsectData#[,-c(218, 222)]                                     

pred_colnames <- c()                                                                       
for (name in names(new_zbjInsectData)) {                                                 
  pred_colnames <- append(pred_colnames, predator_names[[name]])                           
}

colnames(new_zbjInsectData) <- c(pred_colnames)                                          
#####

HR_cols = which(names(new_zbjInsectData) == "Hipposideros ruber")
RA_cols = which(names(new_zbjInsectData) == "Rhinolophus alcyone")

sampleNumbers <- list(HR_nbrs = names(zbjInsectData[,c(HR_cols)]),
                      RA_nbrs = names(zbjInsectData[,c(RA_cols)]))

sampleNumbers <- append(sampleNumbers, 
                        list(allSample_nbrs = c(sampleNumbers$HR_nbrs, sampleNumbers$RA_nbrs)))

focalData <- as.data.frame(zbjInsectData[,c(sampleNumbers$allSample_nbrs)])
rownames(focalData) <- c(insectData$OTU)

# function for converting seq data to occurrence data
##### 
convertToOccurrence <- function(df) { 
  ifelse(df > 0,                    
         df <- 1, df <- 0)
}
#####

# function for converting seq data to RRA data
#####
convertToRRA <- function(df) {
  newDF <- df
  for (row in 1:nrow(df)) {
    for (col in 1:ncol(df)) {
      if (df[row,col] >= 1) {
        newDF[row,col] <- (df[row,col]/sum(df[,col]))*100
      }
    }
  }
  return(newDF)
}
#####

library(dplyr)

focalData_occurrence <- data.frame(lapply(focalData, convertToOccurrence))
colnames(focalData_occurrence) <- c(names(focalData))
rownames(focalData_occurrence) <- c(insectData$OTU)

focalDataList <- list(Original_data = as.data.frame(t(focalData)),
                  Occurrence_data = as.data.frame(t(focalData_occurrence)),
                  RRA_data = as.data.frame(t(convertToRRA(focalData))))

focalInfo <- data.frame()                                  
for (nbr in as.numeric(sampleNumbers$allSample_nbrs)) {                       
  for (row in 1:nrow(zbj_final)) {                          
    ifelse(zbj_final[row, "Lab.nbr"] == nbr,
           focalInfo <- rbind(focalInfo, zbj_final[row,]),
           focalInfo <- focalInfo)                        
  }                                                         
}

focalInfo <- cbind(focalInfo, paste(focalInfo$Location, "-", focalInfo$Season))
colnames(focalInfo) <- c("Lab.nbr", "Animal", "Species", "Location", "Site", "Season", "Location-Season")

focalDataList <- list(Original_data = focalDataList$Original_data <- cbind(focalDataList$Original_data, focalInfo),
                      Occurrence_data = focalDataList$Occurrence_data <- cbind(focalDataList$Occurrence_data, focalInfo),
                      RRA_data = focalDataList$RRA_data <- cbind(focalDataList$RRA_data, focalInfo))

# function for filtering columns
#####
filter_cols <- function(x) {                 
  z <- c()                                          
  for (col in 1:ncol(x)) {                         
    which.col <- which(x[,col] == 0)                
    if (length(which.col)==length(x[,col])) {      
      z <- append(z, col) 
    }                                               
  } 
  return(z)
}
#####

occurrencesCut <- filter_cols(focalDataList$Occurrence_data[1:13975])
RRACut <- filter_cols(focalDataList$RRA_data[1:13975])

focalDataList <- append(focalDataList, 
                        list(Occurrence_data_cut = focalDataList$Occurrence_data[,-c(occurrencesCut)],
                             RRA_data_cut = focalDataList$RRA_data[,-c(RRACut)]))

# total number of reads for all species
#####
readCount <- c()                                            
for (col in 1:ncol(insectData[2:268])) {                  
  readCount <- append(readCount, sum(insectData[2:268][,col]))
}                                                         

allReadCounts <- sum(readCount)
#####

HR_occurrences_data = focalDataList$Occurrence_data_cut[c(which(focalDataList$Occurrence_data_cut$Lab.nbr %in% as.numeric(sampleNumbers$HR_nbrs))),]
RA_occurrences_data = focalDataList$Occurrence_data_cut[c(which(focalDataList$Occurrence_data_cut$Lab.nbr %in% as.numeric(sampleNumbers$RA_nbrs))),]
HR_RRA_data = focalDataList$RRA_data_cut[c(which(focalDataList$RRA_data_cut$Lab.nbr %in% as.numeric(sampleNumbers$HR_nbrs))),]
RA_RRA_data = focalDataList$RRA_data_cut[c(which(focalDataList$RRA_data_cut$Lab.nbr %in% as.numeric(sampleNumbers$RA_nbrs))),]

HR_occurrencesCut <- filter_cols(HR_occurrences_data[1:11342])
RA_occurrencesCut <- filter_cols(RA_occurrences_data[1:11342])
HR_RRACut <- filter_cols(HR_RRA_data[1:11342])
RA_RRACut <- filter_cols(RA_RRA_data[1:11342])

focalDataList <- append(focalDataList, 
                        list(HR_occurrences_data = HR_occurrences_data[,-c(HR_occurrencesCut)],
                             RA_occurrences_data = RA_occurrences_data[,-c(RA_occurrencesCut)],
                             HR_RRA_data = HR_RRA_data[,-c(HR_RRACut)],
                             RA_RRA_data = RA_RRA_data[,-c(RA_RRACut)]))

# taxonomic resolution
#####
focalDataList <- append(focalDataList, list(HR_occurrences_graphs = as.data.frame(t(focalDataList$HR_occurrences_data[1:6501])),
                                            RA_occurrences_graphs = as.data.frame(t(focalDataList$RA_occurrences_data[1:7287])),
                                            HR_RRA_graphs = as.data.frame(t(focalDataList$HR_RRA_data[1:6501])),
                                            RA_RRA_graphs = as.data.frame(t(focalDataList$RA_RRA_data[1:7287]))))

HR_occurrences_OTUs = insectData[c(which(insectData$OTU %in% rownames(focalDataList$HR_occurrences_graphs))),]
RA_occurrences_OTUs = insectData[c(which(insectData$OTU %in% rownames(focalDataList$RA_occurrences_graphs))),]

insectTaxonomy <- list(all_OTUs = insectData[,c(1,269,272:275)],
                       HR_occurrences_OTUs = cbind(focalDataList$HR_occurrences_graphs, HR_occurrences_OTUs[,c(1,269:275)]),
                       RA_occurrences_OTUs = cbind(focalDataList$RA_occurrences_graphs, RA_occurrences_OTUs[,c(1,269:275)]),
                       HR_RRA_OTUs = cbind(focalDataList$HR_RRA_graphs, HR_occurrences_OTUs[,c(1,269:275)]),
                       RA_RRA_OTUs = cbind(focalDataList$RA_RRA_graphs, RA_occurrences_OTUs[,c(1,269:275)]))

not.na <- c()
for (l in 1:length(insectTaxonomy$HR_occurrences_OTUs$family)) {
  ifelse(insectTaxonomy$HR_occurrences_OTUs$family[l] != "NA",
         not.na <- append(not.na, l),
         not.na <- not.na)
}

(length(not.na)/length(insectTaxonomy$HR_occurrences_OTUs$order))*100
#####

# proportion of Lepidoptera
#####
RA.order <- RA.focal.df[-c(which(RA.focal.df$order=="NA")),]
RA_Lepidoptera <- length(which(RA.order$order=="Lepidoptera"))
(RA_Lepidoptera/length(RA.order$order))*100

length(unique(RA.order$family))
#####

# lab numbers for focal species in each location-season
#####
HR_KW <- focalDataList$HR_occurrences_data[focalDataList$HR_occurrences_data$`Location-Season`=="Konye - wet",]
HR_KW_nbrs <- HR_KW$Lab.nbr

HR_KD <- focalDataList$HR_occurrences_data[focalDataList$HR_occurrences_data$`Location-Season`=="Konye - dry",]
HR_KD_nbrs <- HR_KD$Lab.nbr

HR_AW <- focalDataList$HR_occurrences_data[focalDataList$HR_occurrences_data$`Location-Season`=="Ayos - wet",]
HR_AW_nbrs <- HR_AW$Lab.nbr

HR_AD <- focalDataList$HR_occurrences_data[focalDataList$HR_occurrences_data$`Location-Season`=="Ayos - dry",]
HR_AD_nbrs <- HR_AD$Lab.nbr

HR_BW <- focalDataList$HR_occurrences_data[focalDataList$HR_occurrences_data$`Location-Season`=="Bokito - wet",]
HR_BW_nbrs <- HR_BW$Lab.nbr

HR_BD <- focalDataList$HR_occurrences_data[focalDataList$HR_occurrences_data$`Location-Season`=="Bokito - dry",]
HR_BD_nbrs <- HR_BD$Lab.nbr

RA_KW <- focalDataList$RA_occurrences_data[focalDataList$RA_occurrences_data$`Location-Season`=="Konye - wet",]
RA_KW_nbrs <- RA_KW$Lab.nbr

RA_KD <- focalDataList$RA_occurrences_data[focalDataList$RA_occurrences_data$`Location-Season`=="Konye - dry",]
RA_KD_nbrs <- RA_KD$Lab.nbr

RA_AW <- focalDataList$RA_occurrences_data[focalDataList$RA_occurrences_data$`Location-Season`=="Ayos - wet",]
RA_AW_nbrs <- RA_AW$Lab.nbr

RA_AD <- focalDataList$RA_occurrences_data[focalDataList$RA_occurrences_data$`Location-Season`=="Ayos - dry",]
RA_AD_nbrs <- RA_AD$Lab.nbr

RA_BW <- focalDataList$RA_occurrences_data[focalDataList$RA_occurrences_data$`Location-Season`=="Bokito - wet",]
RA_BW_nbrs <- RA_BW$Lab.nbr

RA_BD <- focalDataList$RA_occurrences_data[focalDataList$RA_occurrences_data$`Location-Season`=="Bokito - dry",]
RA_BD_nbrs <- RA_BD$Lab.nbr
#####

# function for filtering rows
#####
filter_rows <- function(x) {                  
     z <- c()                                          
     for (row in 1:nrow(x)) {                          
       which.row <- which(x[row,] == 0)                 
       if (length(which.row)==length(x[row,])) {       
         z <- append(z, row)                           
       }                                               
     }                                                 
     return(z)                                         
   }
#####

#####
HR_KW_occurrence <- insectTaxonomy$HR_occurrences_OTUs[,c(as.character(HR_KW_nbrs))]
HR_KD_occurrence <- insectTaxonomy$HR_occurrences_OTUs[,c(as.character(HR_KD_nbrs))]
HR_AW_occurrence <- insectTaxonomy$HR_occurrences_OTUs[,c(as.character(HR_AW_nbrs))]
HR_AD_occurrence <- insectTaxonomy$HR_occurrences_OTUs[,c(as.character(HR_AD_nbrs))]
HR_BW_occurrence <- insectTaxonomy$HR_occurrences_OTUs[,c(as.character(HR_BW_nbrs))]
HR_BD_occurrence <- insectTaxonomy$HR_occurrences_OTUs[,c(as.character(HR_BD_nbrs))]

RA_KW_occurrence <- insectTaxonomy$RA_occurrences_OTUs[,c(as.character(RA_KW_nbrs))]
RA_KD_occurrence <- insectTaxonomy$RA_occurrences_OTUs[,c(as.character(RA_KD_nbrs))]
RA_AW_occurrence <- insectTaxonomy$RA_occurrences_OTUs[,c(as.character(RA_AW_nbrs))]
RA_AD_occurrence <- insectTaxonomy$RA_occurrences_OTUs[,c(as.character(RA_AD_nbrs))]
RA_BW_occurrence <- insectTaxonomy$RA_occurrences_OTUs[,c(as.character(RA_BW_nbrs))]
RA_BD_occurrence <- insectTaxonomy$RA_occurrences_OTUs[,c(as.character(RA_BD_nbrs))]

HR_KW_RRA <- insectTaxonomy$HR_RRA_OTUs[,c(as.character(HR_KW_nbrs))]
HR_KD_RRA <- insectTaxonomy$HR_RRA_OTUs[,c(as.character(HR_KD_nbrs))]
HR_AW_RRA <- insectTaxonomy$HR_RRA_OTUs[,c(as.character(HR_AW_nbrs))]
HR_AD_RRA <- insectTaxonomy$HR_RRA_OTUs[,c(as.character(HR_AD_nbrs))]
HR_BW_RRA <- insectTaxonomy$HR_RRA_OTUs[,c(as.character(HR_BW_nbrs))]
HR_BD_RRA <- insectTaxonomy$HR_RRA_OTUs[,c(as.character(HR_BD_nbrs))]

RA_KW_RRA <- insectTaxonomy$RA_RRA_OTUs[,c(as.character(RA_KW_nbrs))]
RA_KD_RRA <- insectTaxonomy$RA_RRA_OTUs[,c(as.character(RA_KD_nbrs))]
RA_AW_RRA <- insectTaxonomy$RA_RRA_OTUs[,c(as.character(RA_AW_nbrs))]
RA_AD_RRA <- insectTaxonomy$RA_RRA_OTUs[,c(as.character(RA_AD_nbrs))]
RA_BW_RRA <- insectTaxonomy$RA_RRA_OTUs[,c(as.character(RA_BW_nbrs))]
RA_BD_RRA <- insectTaxonomy$RA_RRA_OTUs[,as.character(RA_BD_nbrs)]

HR_KW_cut <- filter_rows(HR_KW_occurrence)
HR_KD_cut <- filter_rows(HR_KD_occurrence)
HR_AW_cut <- filter_rows(HR_AW_occurrence)
HR_AD_cut <- filter_rows(HR_AD_occurrence)
HR_BW_cut <- filter_rows(HR_BW_occurrence)
HR_BD_cut <- filter_rows(HR_BD_occurrence)

RA_KW_cut <- filter_rows(RA_KW_occurrence)
RA_KD_cut <- filter_rows(RA_KD_occurrence)
RA_AW_cut <- filter_rows(RA_AW_occurrence)
RA_AD_cut <- filter_rows(RA_AD_occurrence)
RA_BW_cut <- which(RA_BW_occurrence[1]==0)
RA_BD_cut <- which(RA_BD_occurrence[1]==0)
#####

insectTaxonomy <- append(insectTaxonomy,
                         list(HR_KW_occurrence = HR_KW_occurrence[-c(HR_KW_cut),],
                              HR_KD_occurrence = HR_KD_occurrence[-c(HR_KD_cut),],
                              HR_AW_occurrence = HR_AW_occurrence[-c(HR_AW_cut),],
                              HR_AD_occurrence = HR_AD_occurrence[-c(HR_AD_cut),],
                              HR_BW_occurrence = HR_BW_occurrence[-c(HR_BW_cut),],
                              HR_BD_occurrence = HR_BD_occurrence[-c(HR_BD_cut),],
                              RA_KW_occurrence = RA_KW_occurrence[-c(RA_KW_cut),],
                              RA_KD_occurrence = RA_KD_occurrence[-c(RA_KD_cut),],
                              RA_AW_occurrence = RA_AW_occurrence[-c(RA_AW_cut),],
                              RA_AD_occurrence = RA_AD_occurrence[-c(RA_AD_cut),],
                              RA_BW_occurrence = RA_BW_occurrence,
                              RA_BD_occurrence = RA_BD_occurrence,
                              HR_KW_RRA = HR_KW_RRA[-c(HR_KW_cut),],
                              HR_KD_RRA = HR_KD_RRA[-c(HR_KD_cut),],
                              HR_AW_RRA = HR_AW_RRA[-c(HR_AW_cut),],
                              HR_AD_RRA = HR_AD_RRA[-c(HR_AD_cut),],
                              HR_BW_RRA = HR_BW_RRA[-c(HR_BW_cut),],
                              HR_BD_RRA = HR_BD_RRA[-c(HR_BD_cut),],
                              RA_KW_RRA = RA_KW_RRA[-c(RA_KW_cut),],
                              RA_KD_RRA = RA_KD_RRA[-c(RA_KD_cut),],
                              RA_AW_RRA = RA_AW_RRA[-c(RA_AW_cut),],
                              RA_AD_RRA = RA_AD_RRA[-c(RA_AD_cut),],
                              RA_BW_RRA = RA_BW_RRA,
                              RA_BD_RRA = RA_BD_RRA,
                              HR_KW_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_KW_cut), c("OTU","sacc","order","family","genus","species")],
                              HR_KD_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_KD_cut), c("OTU","sacc","order","family","genus","species")],
                              HR_AW_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_AW_cut), c("OTU","sacc","order","family","genus","species")],
                              HR_AD_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_AD_cut), c("OTU","sacc","order","family","genus","species")],
                              HR_BW_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_BW_cut), c("OTU","sacc","order","family","genus","species")],
                              HR_BD_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_BD_cut), c("OTU","sacc","order","family","genus","species")],
                              RA_KW_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_KW_cut), c("OTU","sacc","order","family","genus","species")],
                              RA_KD_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_KD_cut), c("OTU","sacc","order","family","genus","species")],
                              RA_AW_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_AW_cut), c("OTU","sacc","order","family","genus","species")],
                              RA_AD_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_AD_cut), c("OTU","sacc","order","family","genus","species")],
                              RA_BW_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_BW_cut), c("OTU","sacc","order","family","genus","species")],
                              RA_BD_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_BD_cut), c("OTU","sacc","order","family","genus","species")]))

# function to count OTUs per sample
#####
countOTUs <- function(df) {                             
  countOTUs <- c()                                   
  for (col in 1:ncol(df)) {                              
    count <- which(df[,col]==1)                          
    countOTUs <- append(countOTUs, length(count))
  }                                                      
  return(countOTUs)                                  
  
} 
#####

# is there a significant difference in the number of OTUs between species
#####
HRcountOTUs <- countOTUs(focalDataList$HR_occurrences_graphs)
RAcountOTUs <- countOTUs(focalDataList$RA_occurrences_graphs)

shapiro.test(c(HRcountOTUs, RAcountOTUs))
wilcox.test(HRcountOTUs, RAcountOTUs)

HR_KW_countOTUs <- countOTUs(insectTaxonomy$HR_KW_occurrence)
HR_KD_countOTUs <- countOTUs(insectTaxonomy$HR_KD_occurrence)
HR_AW_countOTUs <- countOTUs(insectTaxonomy$HR_AW_occurrence)
HR_AD_countOTUs <- countOTUs(insectTaxonomy$HR_AD_occurrence)
HR_BW_countOTUs <- countOTUs(insectTaxonomy$HR_BW_occurrence)
HR_BD_countOTUs <- countOTUs(insectTaxonomy$HR_BD_occurrence)

RA_KW_countOTUs <- countOTUs(insectTaxonomy$RA_KW_occurrence)
RA_KD_countOTUs <- countOTUs(insectTaxonomy$RA_KD_occurrence)
RA_AW_countOTUs <- countOTUs(insectTaxonomy$RA_AW_occurrence)
RA_AD_countOTUs <- countOTUs(insectTaxonomy$RA_AD_occurrence)
RA_BW_countOTUs <- length(which(insectTaxonomy$RA_BW_occurrence ==1))
RA_BD_countOTUs <- length(which(insectTaxonomy$RA_BD_occurrence ==1))

compareOTUs <- data.frame(LS = c(rep("Konye - wet", times = length(RA_KW_countOTUs)),
                                 rep("Konye - dry", times = length(RA_KD_countOTUs)),
                                 rep("Ayos - wet", times = length(RA_AW_countOTUs)),
                                 rep("Ayos - dry", times = length(RA_AD_countOTUs)),
                                 rep("Bokito - wet", times = length(RA_BW_countOTUs)),
                                 rep("Bokito - dry", times = length(RA_BD_countOTUs))),
                          Counts = c(RA_KW_countOTUs, RA_KD_countOTUs,RA_AW_countOTUs,
                                     RA_AD_countOTUs,RA_BW_countOTUs,RA_BD_countOTUs))

compareOTUs$LS <- factor(compareOTUs$LS,
                         levels = c("Konye - wet", "Konye - dry",
                                    "Ayos - wet", "Ayos - dry",
                                    "Bokito - wet", "Bokito - dry"))

library(ggplot2)

compareOTUsBoxplot <- ggplot(compareOTUs, aes(x = LS, y = Counts, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 10),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                                                axis.text.y = element_text(size = 10)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = "Number of OTUs") + ylim(0,300)

compareOTUsBoxplot + annotate("text",
                              x = 1:length(table(compareOTUs$LS)),
                              y = rep(290, times = 6),
                              label = table(compareOTUs$LS),
                              col = "red",
                              vjust = - 1)


compareOTUstats <- data.frame(LS = c(rep("HR - Konye - wet", times = length(HR_KW_countOTUs)),
                                           rep("HR - Konye - dry", times = length(HR_KD_countOTUs)),
                                           rep("HR - Ayos - wet", times = length(HR_AW_countOTUs)),
                                           rep("HR - Ayos - dry", times = length(HR_AD_countOTUs)),
                                           rep("HR - Bokito - wet", times = length(HR_BW_countOTUs)),
                                           rep("HR - Bokito - dry", times = length(HR_BD_countOTUs)),
                                           rep("RA - Konye - wet", times = length(RA_KW_countOTUs)),
                                           rep("RA - Konye - dry", times = length(RA_KD_countOTUs)),
                                           rep("RA - Ayos - wet", times = length(RA_AW_countOTUs)),
                                           rep("RA - Ayos - dry", times = length(RA_AD_countOTUs)),
                                           rep("RA - Bokito - wet", times = length(RA_BW_countOTUs)),
                                           rep("RA - Bokito - dry", times = length(RA_BD_countOTUs))),
                                    Counts = c(HR_KW_countOTUs, HR_KD_countOTUs, HR_AW_countOTUs,
                                               HR_AD_countOTUs, HR_BW_countOTUs, HR_BD_countOTUs,
                                               RA_KW_countOTUs, RA_KD_countOTUs, RA_AW_countOTUs, 
                                               RA_AD_countOTUs, RA_BW_countOTUs, RA_BD_countOTUs))

compareOTUkruskal <- kruskal.test(data = compareOTUstats, x = compareOTUstats$Counts,
                                     g = compareOTUstats$LS, formula = Counts ~ LS)
#####

# occurrence calculations for focal species
#####
wPOOlist <- c()
for (name in names(insectTaxonomy[6:15])) {
  wPOOlist <- append(wPOOlist, wPOO_function(insectTaxonomy[[name]]))
}

for (name in names(insectTaxonomy[30:39])) {
  insectTaxonomy[[name]] <- cbind(insectTaxonomy[[name]], wPOOlist[1:nrow(insectTaxonomy[[name]])])
  wPOOlist <- wPOOlist[-c(1:nrow(insectTaxonomy[[name]]))]
}
#####

# relative read abundance
# function for calculating RRA
#####
RRA_function <- function(df) {
  sumReadCounts <- c()
  
  for (row in 1:nrow(df)) {
    totalReadCounts <- c()
    
    for (col in 1:ncol(df)) {
      
      totalReadCounts <- append(totalReadCounts, as.numeric(df[row,col]/sum(df[,col])))
    }
    sumReadCounts <- append(sumReadCounts, sum(totalReadCounts))
  }
  return(sumReadCounts)
}
#####

#####
RRAlist <- c()
for (name in names(insectTaxonomy[18:27])) {
  sumReadCounts <- RRA_function(insectTaxonomy[[name]])
  RRAlist <- append(RRAlist, sumReadCounts/length(ncol(insectTaxonomy[[name]]))*100)
}

for (name in names(insectTaxonomy[30:39])) {
  insectTaxonomy[[name]] <- cbind(insectTaxonomy[[name]], RRAlist[1:nrow(insectTaxonomy[[name]])])
  RRAlist <- RRAlist[-c(1:nrow(insectTaxonomy[[name]]))]
}
#####

# Shannon diversity 
# function for calculating Shannon index
#####
library(vegan)
ShannonIndexFunction <- function(df) {                                 
  ShannonIndex <- c()                                                   
  
  for (col in 1:ncol(df)) {                                               
    ShannonIndex <- append(ShannonIndex, diversity(df[,col], "shannon"))
  }                                                                       
  return(ShannonIndex)                                                  
} 
#####

#####
HR_occurrences_Shannon <- c()
for (name in names(insectTaxonomy[6:11])) {
  HR_occurrences_Shannon <- append(HR_occurrences_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

HR_RRA_Shannon <- c()
for (name in names(insectTaxonomy[18:23])) {
  HR_RRA_Shannon <- append(HR_RRA_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

RA_occurrences_Shannon <- c()
for (name in names(insectTaxonomy[12:15])) {
  RA_occurrences_Shannon <- append(RA_occurrences_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

RA_RRA_Shannon <- c()
for (name in names(insectTaxonomy[24:27])) {
  RA_RRA_Shannon <- append(RA_RRA_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

compareShannon <- data.frame(LS = c(rep("Konye - wet", times = length(insectTaxonomy$RA_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$RA_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$RA_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$RA_AD_occurrence))),
                                    #rep("Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                    #rep("Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence))),
                                 Shannon.index = RA_occurrences_Shannon)

compareShannon$LS <- factor(compareShannon$LS,
                                levels = c("Konye - wet", "Konye - dry",
                                           "Ayos - wet", "Ayos - dry"))
                                           #"Bokito - wet", "Bokito - dry"))

compareShannonBoxplot <- ggplot(compareShannon, aes(x = LS, y = Shannon.index, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = "Shannon index") + ylim(2.3,7.2)

compareShannonBoxplot + annotate("text",
                                      x = 1:length(table(compareShannon$LS)),
                                      y = rep(6.75, times = 4),
                                      label = table(compareShannon$LS),
                                      col = "red",
                                      vjust = - 1)

compareShannonStats <- data.frame(LS = c(rep("HR - Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                              rep("HR - Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                              rep("HR - Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                              rep("HR - Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                              rep("HR - Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                              rep("HR - Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence)),
                                              rep("RA - Konye - wet", times = length(insectTaxonomy$RA_KW_occurrence)),
                                              rep("RA - Konye - dry", times = length(insectTaxonomy$RA_KD_occurrence)),
                                              rep("RA - Ayos - wet", times = length(insectTaxonomy$RA_AW_occurrence)),
                                              rep("RA - Ayos - dry", times = length(insectTaxonomy$RA_AD_occurrence))),
                                       Counts = c(HR_occurrences_Shannon, RA_occurrences_Shannon))


shapiro.test(compareShannonStats$Counts)
compareShannonkruskal <- kruskal.test(data = compareShannonStats, x = compareShannonStats$Counts,
                                        g = compareShannonStats$LS, formula = Counts ~ LS)

mean(c(HR_occurrences_Shannon, RA_occurrences_Shannon))
wilcox.test(c(HR_occurrences_Shannon, RA_occurrences_Shannon), c(HR_RRA_Shannon, RA_RRA_Shannon))
#####

# Niche breadth
# function for calculating Levin's index
#####
LevinsFunction <- function(df) {
  sumDietProp <- c()
  for (col in 1:ncol(df)) {
    dietProp <- c()
    for (row in 1:nrow(df)) {
      dietProp <- append(dietProp, (df[row,col]/sum(df[,col]))^2)
    }
    sumDietProp <- append(sumDietProp, sum(dietProp))
  }
  
  Levins <- c()
  LevinsStd <- c()
  for (prop in sumDietProp) {
    std <- 1/prop
    Levins <- append(Levins, std)
    LevinsStd <- append(LevinsStd, ((std-1)/(length(df[,1])-1)))
  }
  
  return(LevinsStd)
}
#####

#####
HR_occurrences_Levins <- c()
for (name in names(insectTaxonomy[6:11])) {
  HR_occurrences_Levins <- append(HR_occurrences_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

HR_RRA_Levins <- c()
for (name in names(insectTaxonomy[18:23])) {
  HR_RRA_Levins <- append(HR_RRA_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

RA_occurrences_Levins <- c()
for (name in names(insectTaxonomy[12:15])) {
  RA_occurrences_Levins <- append(RA_occurrences_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

RA_RRA_Levins <- c()
for (name in names(insectTaxonomy[24:27])) {
  RA_RRA_Levins <- append(RA_RRA_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

compareLevins <- data.frame(LS = c(rep("Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                    rep("Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                    rep("Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence))),
                             Levins.index = HR_RRA_Levins)

compareLevins$LS <- factor(compareLevins$LS,
                            levels = c("Konye - wet", "Konye - dry",
                                       "Ayos - wet", "Ayos - dry",
                                       "Bokito - wet", "Bokito - dry"))

y_expression <- expression("Levin's index" ~ "" ~(B[A]))

compareLevinsBoxplot <- ggplot(compareLevins, aes(x = LS, y = Levins.index, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = y_expression) + ylim(0,0.0105)

compareLevinsBoxplot + annotate("text",
                                 x = 1:length(table(compareLevins$LS)),
                                 y = rep(0.0102, times = 6),
                                 label = table(compareLevins$LS),
                                 col = "red",
                                 vjust = - 1)

compareLevinsStats <- data.frame(LS = c(rep("HR - Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                         rep("HR - Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                         rep("HR - Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                         rep("HR - Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                         rep("HR - Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                         rep("HR - Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence)),
                                         rep("RA - Konye - wet", times = length(insectTaxonomy$RA_KW_occurrence)),
                                         rep("RA - Konye - dry", times = length(insectTaxonomy$RA_KD_occurrence)),
                                         rep("RA - Ayos - wet", times = length(insectTaxonomy$RA_AW_occurrence)),
                                         rep("RA - Ayos - dry", times = length(insectTaxonomy$RA_AD_occurrence))),
                                  Counts = c(HR_RRA_Levins, RA_RRA_Levins))


shapiro.test(compareLevinsStats$Counts)
compareShannonkruskal <- kruskal.test(data = compareLevinsStats, x = compareLevinsStats$Counts,
                                      g = compareLevinsStats$LS, formula = Counts ~ LS)

library(FSA)
dunnTest(Counts ~ LS, data = compareLevinsStats, method = "bonferroni")

mean(c(HR_occurrences_Levins, RA_occurrences_Levins))
wilcox.test(c(HR_occurrences_Levins, RA_occurrences_Levins), c(HR_RRA_Levins, RA_RRA_Levins))
#####

# Pianka's niche overlap index
#####
library(EcoSimR)

piankaNullFunction <- function(HRdf, RAdf) {
  nullModel <- c(1:10000)
  
  for (n in nullModel) {
    ra3model <- ra3(speciesData = matrix(c(HRdf, RAdf), nrow = 2))
    
    nullModel[n] <- pianka(m = ra3model) 
  }
  return(nullModel)
}

AD_RRA_pianka <- pianka(m = matrix(c(insectTaxonomy$HR_AD_info$`RRAlist[1:nrow(insectTaxonomy[[name]])]`,
                                         insectTaxonomy$RA_AD_info$`RRAlist[1:nrow(insectTaxonomy[[name]])]`), 
                                       nrow = 2))
AD_RRA_null <- piankaNullFunction(insectTaxonomy$HR_AD_info$`RRAlist[1:nrow(insectTaxonomy[[name]])]`,
                                      insectTaxonomy$RA_AD_info$`RRAlist[1:nrow(insectTaxonomy[[name]])]`)
wilcox.test(AD_RRA_null, mu = AD_RRA_pianka)

LocationSeasonList <- c("Konye - wet" , "Konye - dry" , "Ayos - wet" , "Ayos - dry")

piankaData <- data.frame(LS = rep(LocationSeasonList, times = 2),
                          Metric = c(rep("wPOO", times = 4), rep("RRA", times = 4)),
                          Pianka = c(KW_wPOO_pianka, KD_wPOO_pianka,
                                     AW_wPOO_pianka, AD_wPOO_pianka,
                                     KW_RRA_pianka, KD_RRA_pianka,
                                     AW_RRA_pianka, AD_RRA_pianka)) 

piankaData$LS <- factor(piankaData$LS,
                         levels = c("Konye - wet", "Konye - dry",
                                    "Ayos - wet", "Ayos - dry"))

piankaPlot <- ggplot(piankaData, aes(x = LS, y = Pianka, fill = Metric)) + 
  geom_bar(stat="identity", position = "dodge") +
  theme_minimal() + scale_fill_manual(values = c("darkolivegreen3", "darkorange")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  labs(x = "Location - Season", y = "Pianka's index")

hist(piankaData$Pianka)
shapiro.test(piankaData$Pianka)
t.test(c(KW_wPOO_pianka, KD_wPOO_pianka,
         AW_wPOO_pianka, AD_wPOO_pianka), 
       c(KW_RRA_pianka, KD_RRA_pianka,
         AW_RRA_pianka, AD_RRA_pianka))

wilcox.test(c(HR.KW.wPOO.pianka, HR.KD.wPOO.pianka,
              HR.AW.wPOO.pianka, HR.AD.wPOO.pianka), 
            c(HR.KW.RRA.pianka, HR.KD.RRA.pianka,
              HR.AW.RRA.pianka, HR.AD.RRA.pianka))

wPOO.Pianka <- pianka.data[pianka.data$Metric=="RRA",]
mean(c(KW_RRA_pianka, KD_RRA_pianka,
       AW_RRA_pianka, AD_RRA_pianka))
#####

# db-rda
# with presence/absence data
capscaleFunction <- function(df, distance) {
  
  focalInteract <- capscale(df[1:11336] ~ df$Species, distance = distance, add = TRUE)
  interactOutput <- anova(focalInteract, step = 1000, perm.max = 200)
  
  vec = interactOutput$SumOfSqs/sum(interactOutput$SumOfSqs)*100
  table = interactOutput
  table$SumOfSqs = vec
  table
  
  return(table)
}

LocationSeasonList <- c("Konye - wet", "Konye - dry", "Ayos - wet",
                        "Ayos - dry", "Bokito - wet", "Bokito - dry")


for (LS in LocationSeasonList) {
  print(capscaleFunction(focalDataList$Occurrence_data[focalDataList$Occurrence_data$`Location-Season`==LS,], "jaccard"))
  
  print(capscaleFunction(focalDataList$RRA_data[focalDataList$RRA_data$`Location-Season`==LS,], "bray"))
}

# insect orders present in the diet
unique(insectTaxonomy$HR_occurrences_OTUs$order)

insectTaxonomy[["HR_occurrences_OTUs"]] <- cbind(insectTaxonomy[["HR_occurrences_OTUs"]],
                                                 wPOO_function(insectTaxonomy[["HR_occurrences_OTUs"]][1:77]))

insectTaxonomy[["RA_occurrences_OTUs"]] <- cbind(insectTaxonomy[["RA_occurrences_OTUs"]],
                                                 wPOO_function(insectTaxonomy[["RA_occurrences_OTUs"]][1:82]))

# function to find most common Lepidopteran families in each location and season
findLepidoptera <- function(df) {
  lepidoptera_df <- df[df$order=="Lepidoptera",]
  lepidoptera_fam <- unique(lepidoptera_df$family)
  fam_na <- which(lepidoptera_fam=="NA")
  
  return(lepidoptera_fam[-fam_na])
}

HRorder <- insectTaxonomy$HR_occurrences_OTUs[-c(which(insectTaxonomy$HR_occurrences_OTUs$order=="NA")),]
HRpropLepidoptera <- length(which(HRorder$order=="Lepidoptera"))/length(HRorder$order)

HRLepFam <- length(unique(HRorder$family))

RAorder <- insectTaxonomy$RA_occurrences_OTUs[-c(which(insectTaxonomy$RA_occurrences_OTUs$order=="NA")),]
RApropLepidoptera <- length(which(RAorder$order=="Lepidoptera"))/length(RAorder$order)

RALepFam <- length(unique(RAorder$family))

# find Lepidopteran families in all location-seasons for each species
LepFam <- list()
for (name in names(insectTaxonomy[30:41])) {
  LepFam <- append(LepFam, list(findLepidoptera(insectTaxonomy[[name]])))
}

names(LepFam) <- c("HR_KW_fam", "HR_KD_fam", "HR_AW_fam",
                   "HR_AD_fam", "HR_BW_fam", "HR_BD_fam",
                   "RA_KW_fam", "RA_KD_fam", "RA_AW_fam", 
                   "RA_AD_fam", "RA_BW_fam", "RA_BD_fam")

allLepidoptera <- findLepidoptera(insectTaxonomy$all_OTUs)

lepidopteraMatrix <- matrix(ncol = 0, nrow = length(allLepidoptera))

for (name in names(LepFam[1:6])) {
  yesNo <- c()
  
  for (all in allLepidoptera) {
    ifelse(all %in% LepFam[[name]], 
           yesNo <- append(yesNo, 1,),
           yesNo <- append(yesNo, 0))
  }
  lepidopteraMatrix <- cbind(lepidopteraMatrix, yesNo)
} 

for (name in names(LepFam[7:12])) {
  yesNo <- c()
  
  for (all in allLepidoptera) {
    ifelse(all %in% LepFam[[name]], 
           yesNo <- append(yesNo, 2,),
           yesNo <- append(yesNo, 0))
  }
  lepidopteraMatrix <- cbind(lepidopteraMatrix, yesNo)
}


colnames(lepidopteraMatrix) <- c(LocationSeasonList)
rownames(lepidopteraMatrix) <- c(allLepidoptera)

library("pheatmap")
pheatmap(lepidopteraMatrix, border_color = "white",
         cluster_rows = FALSE, cluster_cols = FALSE,  angle_col = 0,
         legend = TRUE)
