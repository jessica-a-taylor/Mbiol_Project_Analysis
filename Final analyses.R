# import seq data
seqData <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_consensus (1).xlsx")

# filter for insect OTUs
seqInsecta <- which(seqData$class=="Insecta")
insectData <- seq.data[c(seqInsecta),]

# remove unwanted columns
insectData <- insectData[,-c(56,185,248,249,272,274:275,277:291)]
insectData <- cbind(insectData, paste(insectData$order, insectData$sacc, sep = ""))

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
                       HR_occurrences_OTUs = cbind(focalDataList$HR_occurrences_graphs, HR_occurrences_OTUs[,c(1,269:275, 277)]),
                       RA_occurrences_OTUs = cbind(focalDataList$RA_occurrences_graphs, RA_occurrences_OTUs[,c(1,269:275, 277)]),
                       HR_RRA_OTUs = cbind(focalDataList$HR_RRA_graphs, HR_occurrences_OTUs[,c(1,269:275, 277)]),
                       RA_RRA_OTUs = cbind(focalDataList$RA_RRA_graphs, RA_occurrences_OTUs[,c(1,269:275, 277)]))

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
                              HR_KW_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_KW_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              HR_KD_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_KD_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              HR_AW_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_AW_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              HR_AD_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_AD_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              HR_BW_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_BW_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              HR_BD_info = insectTaxonomy$HR_occurrences_OTUs[-c(HR_BD_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              RA_KW_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_KW_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              RA_KD_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_KD_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              RA_AW_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_AW_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              RA_AD_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_AD_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              RA_BW_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_BW_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")],
                              RA_BD_info = insectTaxonomy$RA_occurrences_OTUs[-c(RA_BD_cut), c("OTU","sacc","order","family","genus","species", "paste(insectData$order, insectData$sacc, sep = \"\")")]))

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
  RRA <- c()
  
  for (row in 1:nrow(df)) {
    totalReadCounts <- c()
    
    for (col in 1:ncol(df)) {
      
      totalReadCounts <- append(totalReadCounts, as.numeric(df[row,col]/sum(df[,col])))
    }
    RRA <- append(RRA, (1/ncol(df))*sum(totalReadCounts)*100)
  }
  return(RRA)
}
#####

#####
RRAlist <- c()
for (name in names(insectTaxonomy[18:27])) {
  RRAlist <- append(RRAlist,RRA_function(insectTaxonomy[[name]]))
}

for (name in names(insectTaxonomy[30:39])) {
  insectTaxonomy[[name]] <- cbind(insectTaxonomy[[name]], RRAlist[c(1:nrow(insectTaxonomy[[name]]))])
  RRAlist <- RRAlist[-c(1:nrow(insectTaxonomy[[name]]))]
}
#####

# summarise each occurrence as wPOO or RRA
wPOO_list <- list()
RRA_list <- list()
for (name in names(insectTaxonomy[30:39])) {
  wPOO_list <- append(wPOO_list, list(insectTaxonomy[[name]]$`wPOOlist[1:nrow(insectTaxonomy[[name]])]`))
  RRA_list <- append(RRA_list, list(insectTaxonomy[[name]]$`RRAlist[c(1:nrow(insectTaxonomy[[name]]))]`))
}

names(wPOO_list) <- c(names(insectTaxonomy[6:15]))
names(RRA_list) <- c(names(insectTaxonomy[6:15]))

for (name in names(insectTaxonomy[6:15])) {
  df <- insectTaxonomy[[name]]
  for (row in 1:nrow(df)) {
    for (col in 1:ncol(df)) {
      ifelse(df[row,col]==1,
             df[row,col] <- wPOO_list[[name]][row],
             df[row,col] <- df[row,col])
    }
  }
  insectTaxonomy <- append(insectTaxonomy, list(df))
}

names(insectTaxonomy) <- c(names(insectTaxonomy[1:41]), "HR_KW_wPOOsummary",
                           "HR_KD_wPOOsummary", "HR_AW_wPOOsummary", 
                           "HR_AD_wPOOsummary", "HR_BW_wPOOsummary",
                           "HR_BD_wPOOsummary", "RA_KW_wPOOsummary",
                           "RA_KD_wPOOsummary", "RA_AW_wPOOsummary",
                           "RA_AD_wPOOsummary")

for (name in names(insectTaxonomy[6:15])) {
  df <- insectTaxonomy[[name]]
  for (row in 1:nrow(df)) {
    for (col in 1:ncol(df)) {
      ifelse(df[row,col]==1,
             df[row,col] <- RRA_list[[name]][row],
             df[row,col] <- df[row,col])
    }
  }
  insectTaxonomy <- append(insectTaxonomy, list(df))
}

names(insectTaxonomy) <- c(names(insectTaxonomy[1:51]), "HR_KW_RRAsummary",
                           "HR_KD_RRAsummary", "HR_AW_RRAsummary", 
                           "HR_AD_RRAsummary", "HR_BW_RRAsummary",
                           "HR_BD_RRAsummary", "RA_KW_RRAsummary",
                           "RA_KD_RRAsummary", "RA_AW_RRAsummary",
                           "RA_AD_RRAsummary")

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
for (name in names(insectTaxonomy[42:47])) {
  HR_occurrences_Shannon <- append(HR_occurrences_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

HR_RRA_Shannon <- c()
for (name in names(insectTaxonomy[52:57])) {
  HR_RRA_Shannon <- append(HR_RRA_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

RA_occurrences_Shannon <- c()
for (name in names(insectTaxonomy[48:51])) {
  RA_occurrences_Shannon <- append(RA_occurrences_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

RA_RRA_Shannon <- c()
for (name in names(insectTaxonomy[58:61])) {
  RA_RRA_Shannon <- append(RA_RRA_Shannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

wPOOonly_Shannon <- c()
for (name in names(insectTaxonomy[30:39])) {
  wPOOonly_Shannon <- append(wPOOonly_Shannon, diversity(insectTaxonomy[[name]]$wPOOlist[1:nrow(insectTaxonomy[[name]])], "shannon"))
}


RRAonly_Shannon <- c()
for (name in names(insectTaxonomy[30:39])) {
  RRAonly_Shannon <- append(RRAonly_Shannon, diversity(insectTaxonomy[[name]]$`RRAlist[c(1:nrow(insectTaxonomy[[name]]))]`, "shannon"))
}


compareShannon <- data.frame(Species = c(rep("H.ruber", times = length()),
                                         rep("R.alcyone", times = length())),
                             LS = c(rep("Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                    rep("Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                    rep("Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence)),
                                    rep("Konye - wet", times = length(insectTaxonomy$RA_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$RA_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$RA_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$RA_AD_occurrence))),
                                 Shannon.index = c(HR_occurrences_Shannon,RA_occurrences_Shannon))

compareShannon$LS <- factor(compareShannon$LS,
                                levels = c("Konye - wet", "Konye - dry",
                                           "Ayos - wet", "Ayos - dry",
                                           "Bokito - wet", "Bokito - dry"))

compareShannon <- data.frame(Species = c(rep("H.ruber", times = length(HR_occurrences_Shannon)),
                                         rep("R.alcyone", times = length(RA_occurrences_Shannon))),
                             LS = c(rep("Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                    rep("Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                    rep("Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence)),
                                    rep("Konye - wet", times = length(insectTaxonomy$RA_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$RA_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$RA_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$RA_AD_occurrence))),
                             Shannon.index = c(HR_occurrences_Shannon,RA_occurrences_Shannon))


compareShannonBoxplot <- ggplot(compareShannon,aes(x=LS, y=Shannon.index, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + ylim(2,7) +
  labs(x = "Location - Season", y = "Shannon index") + facet_wrap(vars(Species)) +
  geom_point(aes(x = LS, y = wPOO))

compareShannonPoint <- ggplot(compareShannon, aes(x = LS, y = wPOO)) + 
  geom_point()

compareShannonBoxplot + annotate("text",
                                      x = 1:length(table(compareShannon$LS)),
                                      y = rep(3, times = 4),
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

library(FSA)
dunnTest(Counts ~ LS, data = compareShannonStats, method = "bonferroni")

mean(c(HR_occurrences_Shannon, RA_occurrences_Shannon))
wilcox.test(c(HR_occurrences_Shannon, RA_occurrences_Shannon), c(HR_RRA_Shannon, RA_RRA_Shannon))
#####

# Niche breadth
# function for calculating Levin's index on incidence data
#####
LevinsFunction <- function(df) {
  Levins <- c()
  std <- c()
  for (col in 1:ncol(df)) {
    dietProp <- c()
    for (row in 1:nrow(df)) {
             dietProp <- append(dietProp, (df[row,col]/sum(df[,col]))^2)
    }
    Levins <- append(Levins, 1/sum(dietProp))
  
    std <- append(std, length(which(df[,col]>0))-1)
  }
  
  LevinsStd <- (Levins-1)/std
  
  
  return(LevinsStd)
}
#####
HR_occurrences_Levins <- c()
for (name in names(insectTaxonomy[42:47])) {
  HR_occurrences_Levins <- append(HR_occurrences_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

HR_RRA_Levins <- c()
for (name in names(insectTaxonomy[52:57])) {
  HR_RRA_Levins <- append(HR_RRA_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

RA_occurrences_Levins <- c()
for (name in names(insectTaxonomy[48:51])) {
  RA_occurrences_Levins <- append(RA_occurrences_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

RA_RRA_Levins <- c()
for (name in names(insectTaxonomy[58:61])) {
  RA_RRA_Levins <- append(RA_RRA_Levins, LevinsFunction(insectTaxonomy[[name]]))
}

compareLevins <- data.frame(LS = c(rep("Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                    rep("Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                    rep("Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                    rep("Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                    rep("Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                    rep("Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence))),
                             Levins.index = HR_occurrences_Levins)

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

compareLevinsStats <- data.frame(Species = c(rep("H.ruber", times = length(HR_occurrences_Levins)),
                                             rep("R.alcyone", times = length(RA_occurrences_Levins))),
                                 LS = c(rep("HR - Konye - wet", times = length(insectTaxonomy$HR_KW_occurrence)),
                                         rep("HR - Konye - dry", times = length(insectTaxonomy$HR_KD_occurrence)),
                                         rep("HR - Ayos - wet", times = length(insectTaxonomy$HR_AW_occurrence)),
                                         rep("HR - Ayos - dry", times = length(insectTaxonomy$HR_AD_occurrence)),
                                         rep("HR - Bokito - wet", times = length(insectTaxonomy$HR_BW_occurrence)),
                                         rep("HR - Bokito - dry", times = length(insectTaxonomy$HR_BD_occurrence)),
                                         rep("RA - Konye - wet", times = length(insectTaxonomy$RA_KW_occurrence)),
                                         rep("RA - Konye - dry", times = length(insectTaxonomy$RA_KD_occurrence)),
                                         rep("RA - Ayos - wet", times = length(insectTaxonomy$RA_AW_occurrence)),
                                         rep("RA - Ayos - dry", times = length(insectTaxonomy$RA_AD_occurrence))),
                                  Counts = c(HR_occurrences_Levins, RA_occurrences_Levins))


shapiro.test(compareLevinsStats$Counts)
compareLevinsruskal <- kruskal.test(data = compareLevinsStats, 
                                    x = compareLevinsStats$Counts,
                                    g = compareLevinsStats$LS, 
                                    formula = Counts ~ LS)

library(FSA)
dunnTest(Counts ~ LS, data = compareLevinsStats, method = "bonferroni")

mean(c(HR_RRA_Levins, RA_RRA_Levins[-which(RA_RRA_Levins=="NaN")]))
wilcox.test(c(HR_occurrences_Levins, RA_occurrences_Levins[-which(RA_occurrences_Levins=="NaN")]), c(HR_RRA_Levins, RA_RRA_Levins[-which(RA_RRA_Levins=="NaN")]))
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

lepidopteraMatrix <- lepidopteraMatrix[,1:6] + lepidopteraMatrix[,7:12]
LocationSeasonList <- c("Konye - wet" , "Konye - dry" , "Ayos - wet" , "Ayos - dry", "Bokito - wet", "Bokito - dry")

colnames(lepidopteraMatrix) <- c(LocationSeasonList)
rownames(lepidopteraMatrix) <- c(allLepidoptera)

library("pheatmap")
pheatmap(lepidopteraMatrix, border_color = "white", 
         color = c("white", "khaki2", "indianred", "darkseagreen"),
         cluster_rows = FALSE, cluster_cols = FALSE,  angle_col = 0,
         legend = TRUE)

famSum <- c()
for (row in 1:nrow(lepidopteraMatrix)) {
  famSum <- append(famSum, sum(lepidopteraMatrix[row,]))
}

famSumDF <- data.frame(Family = allLepidoptera,
                       Sum = famSum)
famSumDF <- famSumDF[order(famSumDF[,2]),]

allLepidoptera <- famSumDF$Family
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

lepidopteraMatrix <- lepidopteraMatrix[,1:6] + lepidopteraMatrix[,7:12]

colnames(lepidopteraMatrix) <- c(LocationSeasonList)
rownames(lepidopteraMatrix) <- c(allLepidoptera)
lepidopteraMatrix <- lepidopteraMatrix[-c(1:2),]

pheatmap(lepidopteraMatrix, border_color = "white", 
         color = c("peachpuff1", "lightpink3", "indianred", "coral4"),
         cluster_rows = FALSE, cluster_cols = FALSE,  angle_col = 0,
         legend = TRUE, legend_labels = c("None", "HR only", "RA only", "Both"),
         legend_breaks = c(0,1,2,3))


# pest occurrences
pestData <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\Pest sequences.xlsx")
pestData <- as.data.frame(pestData)
pestData <- pestData[,-c(56,185,272,274:275, 277:291)]

rownames(pestData) <- c(pestData$OTU)

pestData <- cbind(pestData, paste(pestData$Pest,pestData$sacc, sep = ""))

pestDataList <- list()
pestNbrs <- list()
for (name in names(insectTaxonomy[30:41])) {
  pestDataList <- append(pestDataList, list(insectTaxonomy[[name]][c(which(rownames(insectTaxonomy[[name]]) %in% pestData$OTU)),]))
  pestNbrs <- append(pestNbrs, list(which(rownames(insectTaxonomy[[name]]) %in% pestData$OTU)))
}

names(pestNbrs) <- c(names(insectTaxonomy[6:17]))

for (name in names(insectTaxonomy[6:15])) {
  pestDataList <- append(pestDataList, list(insectTaxonomy[[name]][c(pestNbrs[[name]]),]))
}

for (name in names(insectTaxonomy[16:17])) {
  pestDataList <- append(pestDataList, list(insectTaxonomy[[name]][c(pestNbrs[[name]])]))
}

names(pestDataList) <- c("HR_KW_info", "HR_KD_info", "HR_AW_info",
                         "HR_AD_info", "HR_BW_info", "HR_BD_info",
                         "RA_KW_info", "RA_KD_info", "RA_AW_info", 
                         "RA_AD_info", "RA_BW_info", "RA_BD_info",
                         "HR_KW_pest", "HR_KD_pest", "HR_AW_pest",
                         "HR_AD_pest", "HR_BW_pest", "HR_BD_pest",
                         "RA_KW_pest", "RA_KD_pest", "RA_AW_pest", 
                         "RA_AD_pest", "RA_BW_pest", "RA_BD_pest")

pestCut <- list()
for (name in names(pestDataList[23:24])) {
  pestCut <- append(pestCut, list(which(pestDataList[[name]]==0)))
}

pestDataList[["RA_BW_info"]] <- pestDataList[["RA_BW_info"]][-c(pestCut[[1]]),]
pestDataList[["RA_BW_pest"]] <- pestDataList[["RA_BW_pest"]][-c(pestCut[[1]])]
pestDataList[["RA_BD_info"]] <- pestDataList[["RA_BD_info"]][-c(pestCut[[2]]),]
pestDataList[["RA_BD_pest"]] <- pestDataList[["RA_BD_pest"]][-c(pestCut[[2]])]

allPests <- unique(c(pestData$`paste(pestData$Pest, pestData$sacc, sep = "")`))

findPests <- function(df) {
  pestDF <- pestData[c(which(pestData$OTU %in% df$OTU)),]
  pestSacc <- pestDF$`paste(pestData$Pest, pestData$sacc, sep = "")`
  return(pestSacc)
}

for (name in names(pestDataList[1:10])) {
  pestDataList[[name]] <- cbind(pestDataList[[name]], findPests(pestDataList[[name]]))
}

names(pestDataList) <- c("HR_KW_info", "HR_KD_info", "HR_AW_info",
                         "HR_AD_info", "HR_BW_info", "HR_BD_info",
                         "RA_KW_info", "RA_KD_info", "RA_AW_info", 
                         "RA_AD_info", "RA_BW_info", "RA_BD_info",
                         "HR_KW_pest", "HR_KD_pest", "HR_AW_pest",
                         "HR_AD_pest", "HR_BW_pest", "HR_BD_pest",
                         "RA_KW_pest", "RA_KD_pest", "RA_AW_pest", 
                         "RA_AD_pest", "RA_BW_pest", "RA_BD_pest")

# matrix of pest occurrences
pestMatrix <- matrix(ncol = 0, nrow = length(allPests))

for (name in names(pestDataList[7:10])) {
  yesNo <- c()
  
  for (all in allPests) {
    ifelse(all %in% pestDataList[[name]]$`findPests(pestDataList[[name]])`, 
           yesNo <- append(yesNo, 1),
           yesNo <- append(yesNo, 0))
  }
  pestMatrix <- cbind(pestMatrix, yesNo)
} 

rownames(pestMatrix) <- c(allPests)
colnames(pestMatrix) <- c(names(pestDataList[1:4]))


pestMatrixCut <- filter_rows(pestMatrix)
pestMatrix <- pestMatrix[-c(pestMatrixCut),]

rownames(pestMatrix) <- c(abbreviate(c(rownames(pestMatrix)), minlength = 3, 
                                     strict = TRUE, method = "left.keep", use.classes = FALSE))

# pest networks
pestAbbr <- c("Lep", "Col", "Dip", "Hem")

pestTotal <- list()
for (pest in pestAbbr) {
  pestDF <- data.frame(matrix(pestMatrix[which(str_detect(pest,rownames(pestMatrix))==TRUE),], ncol = 4))
  
  pestList <- c()
  for (col in 1:ncol(pestDF)) {
    pestList <- append(pestList, sum(pestDF[,col])) 
  }
  pestTotal <- append(pestTotal, list(pestList))
}

library(bipartite)

pestBipartite <- matrix(c(pestTotal[[1]], pestTotal[[2]], pestTotal[[3]],pestTotal[[4]]), ncol=4)
rownames(pestBipartite) <- c(LocationSeasonList[-c(5:6)])
colnames(pestBipartite) <- c(pestAbbr)

plotweb(pestBipartite, 
        text.rot=90, labsize = 1.2,  
        y.width.low=0.05,y.width.high=0.05, 
        col.high="light blue",col.low="light green",
        col.interaction = "grey90", ybig = 1.5,
        x.lim = c(0,1.5), sequence = list(seq.low = c(LocationSeasonList[-c(5:6)]),
                                          seq.high = c(pestAbbr)), y.lim = c(-0.5,2)) 


visweb(pestBipartite)


# Venn diagram comparing pest consumption
library(RAM)
  # comparing total pest consumption between species
group.venn(list(H.ruber = unique(c(pestDataList$HR_KW_info$`findPests(pestDataList[[name]])`,
                                     pestDataList$HR_KD_info$`findPests(pestDataList[[name]])`,
                                     pestDataList$HR_AW_info$`findPests(pestDataList[[name]])`,
                                     pestDataList$HR_AD_info$`findPests(pestDataList[[name]])`,
                                   pestDataList$HR_BW_info$`findPests(pestDataList[[name]])`,
                                   pestDataList$HR_BD_info$`findPests(pestDataList[[name]])`)),
                    R.alcyone = unique(c(pestDataList$RA_KW_info$`findPests(pestDataList[[name]])`,
                                       pestDataList$RA_KD_info$`findPests(pestDataList[[name]])`,
                                       pestDataList$RA_AW_info$`findPests(pestDataList[[name]])`,
                                       pestDataList$RA_AD_info$`findPests(pestDataList[[name]])`,
                                       pestDataList$RA_BW_info$`findPests(pestDataList[[name]])`,
                                       pestDataList$RA_BD_info$`findPests(pestDataList[[name]])`))),
           label = FALSE, lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(330, 380), cat.dist = 0.05)

# make better Venn function
#####
JessVenn <- function (vectors, cat.cex = 1.5, cex = 1, cat.pos = NULL, cat.dist = NULL, 
                      label = TRUE, lab.cex = 1, lab.col = "black", fill = NULL, 
                      file = NULL, ext = NULL, width = 8, height = 8) 
{
  save <- !is.null(file)
  if (save) {
    .get.dev(file, ext, height = height, width = width)
  }
  if (!requireNamespace("VennDiagram")) {
    stop("package 'VennDiagram' is required for this function")
  }
  if (!requireNamespace("RColorBrewer")) {
    stop("package 'RColorBrewer' is required for this function")
  }
  if (!requireNamespace("grid")) {
    stop("package 'grid' is required to use this function")
  }
  len <- length(vectors)
  if (is.null(fill)) {
    if (len == 2) {
      fill = c("lightpink", "lightblue")
    }
    else {
      fill = RColorBrewer::brewer.pal(len, "Pastel1")
    }
  }
  else {
    if (length(fill) == len) {
      fill = fill
    }
    else if (length(fill) > len) {
      warning(paste("more colors being provided than required, will ignore ", 
                    length(fill) - len, " colors", sep = ""))
      fill = fill[1:len]
    }
    else {
      warning("not enough colors being provided, will use default")
      if (len == 2) {
        fill = c("lightpink", "lightblue")
      }
      else {
        fill = RColorBrewer::brewer.pal(len, "Pastel1")
      }
    }
  }
  if (len > 2 && label) {
    warning("currently only support 2 groups to have actual item labels; will only use numbers")
  }
  else if (len > 5 || len < 2) {
    stop("please provide 2 to 5 vectors")
  }
  alpha = rep(0.5, len)
  if (!is.null(cat.pos) && !is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.dist = cat.dist, cat.pos = cat.pos, 
                                   cat.fontface = "bold", cat.cex = cat.cex, cex = cex, 
                                   filename = NULL, inverted = TRUE)
  }
  else if (!is.null(cat.pos) && is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.pos = cat.pos, cat.fontface = "bold", 
                                   cat.cex = cat.cex, cex = cex, filename = NULL)
  }
  else if (is.null(cat.pos) && !is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.fontface = "bold", cat.dist = cat.dist, 
                                   cat.cex = cat.cex, cex = cex, filename = NULL)
  }
  else {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.fontface = "bold", cat.cex = cat.cex, 
                                   cex = cex, filename = NULL)
  }
  if (len > 2 && len <= 5) {
    grid::grid.newpage()
    grid::grid.draw(v)
  }
  if (len == 2) {
    if (!label) {
      grid::grid.newpage()
      grid::grid.draw(v)
    }
    else {
      name <- lapply(v, names)
      lapply(v, function(i) i$label)
      v.labels <- lapply(v, function(i) i$label)
      v.lab <- vector()
      for (i in 1:length(v.labels)) {
        if (length(v.labels[[i]] %in% names(vectors)) != 
            0 && isTRUE(v.labels[[i]] %in% names(vectors))) {
          v.lab <- c(v.lab, v.labels[[i]])
        }
      }
      v1 <- vectors[[v.lab[1]]]
      v2 <- vectors[[v.lab[2]]]
      v[[5]]$label <- paste(c(v[[5]]$label, setdiff(v1, 
                                                    v2)), collapse = "\n")
      v[[5]]$gp$cex <- lab.cex
      v[[5]]$gp$col <- lab.col
      v[[6]]$label <- paste(c(v[[6]]$label, setdiff(v2, 
                                                    v1)), collapse = "\n")
      v[[6]]$gp$cex <- lab.cex
      v[[6]]$gp$col <- lab.col
      v[[7]]$label <- paste(c(v[[7]]$label, intersect(v1, 
                                                      v2)), collapse = "\n")
      v[[7]]$gp$cex <- lab.cex
      v[[7]]$gp$col <- lab.col
      grid::grid.newpage()
      grid::grid.draw(v)
    }
  }
  if (save) {
    dev.off()
  }
  invisible()
}
#####
  # comparing pest consumption by each species between seasons
JessVenn(list(Wet = unique(c(pestDataList$RA_KW_info$`findPests(pestDataList[[name]])`,
                               pestDataList$RA_AW_info$`findPests(pestDataList[[name]])`)),
           Dry = unique(c(pestDataList$RA_KD_info$`findPests(pestDataList[[name]])`,
                                     pestDataList$RA_AD_info$`findPests(pestDataList[[name]])`))),
                label = FALSE, lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(220, 150), cat.dist = 0.05, fill = c("lightblue2", "lightpink"))

group.venn(list(KonyeWet = unique(pestDataList$HR_KW_info$`findPests(pestDataList[[name]])`),
                KonyeDry = unique(pestDataList$HR_AW_info$`findPests(pestDataList[[name]])`),
                AyosWet = unique(pestDataList$HR_KD_info$`findPests(pestDataList[[name]])`),
                AyosDry = unique(pestDataList$HR_AD_info$`findPests(pestDataList[[name]])`)),
           label = FALSE, lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.dist = c(0.23,0.23,0.12,0.12))

JessVenn(list(HRwet = unique(c(pestDataList$HR_KW_info$`findPests(pestDataList[[name]])`,
                                 pestDataList$HR_AW_info$`findPests(pestDataList[[name]])`)),
                HRdry = unique(c(pestDataList$HR_KD_info$`findPests(pestDataList[[name]])`,
                                 pestDataList$HR_AD_info$`findPests(pestDataList[[name]])`)),
                RAwet = unique(c(pestDataList$RA_KW_info$`findPests(pestDataList[[name]])`,
                                 pestDataList$RA_AW_info$`findPests(pestDataList[[name]])`)),
                RAdry = unique(c(pestDataList$RA_KD_info$`findPests(pestDataList[[name]])`,
                                 pestDataList$RA_AD_info$`findPests(pestDataList[[name]])`))),
           label = FALSE, lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.dist = c(0.23,0.23,0.12,0.12))
#####
RApestMatrix <- matrix(ncol = 0, nrow = length(allPests))

for (name in names(pestDataList[7:12])) {
  yesNo <- c()
  
  for (all in allPests) {
    ifelse(all %in% pestDataList[[name]]$`findPests(pestDataList[[name]])`, 
           yesNo <- append(yesNo, 2),
           yesNo <- append(yesNo, 0))
  }
  RApestMatrix <- cbind(RApestMatrix, yesNo)
} 


comparePestMatrix <- HRpestMatrix[,1:6] + RApestMatrix[,1:6]
rownames(comparePestMatrix) <- c(allPests)
colnames(comparePestMatrix) <- c(LocationSeasonList)

pestSum <- c()
for (row in 1:nrow(comparePestMatrix)) {
  pestSum <- append(pestSum, sum(comparePestMatrix[row,]))
}

pestSumDF <- data.frame(Family = allPests,
                        Sum = pestSum)
pestSumDF <- pestSumDF[order(pestSumDF[,2]),]

allPests <- pestSumDF$Family

RAcomparePestMatrix <- matrix(ncol = 0, nrow = length(allPests))

for (name in names(pestDataList[7:12])) {
  yesNo <- c()
  
  for (all in allPests) {
    ifelse(all %in% pestDataList[[name]]$`findPests(pestDataList[[name]])`, 
           yesNo <- append(yesNo, 2),
           yesNo <- append(yesNo, 0))
  }
  RAcomparePestMatrix <- cbind(RAcomparePestMatrix, yesNo)
} 


comparePestMatrix <- HRcomparePestMatrix[,1:6] + RAcomparePestMatrix[,1:6]
rownames(comparePestMatrix) <- c(allPests)
colnames(comparePestMatrix) <- c(LocationSeasonList)

pestMatrixCut <- filter_rows(comparePestMatrix)
comparePestMatrix <- comparePestMatrix[-c(pestMatrixCut),]

pheatmap(comparePestMatrix, border_color = "white", 
         color = c("peachpuff1", "lightpink3", "indianred", "coral4"),
         cluster_rows = FALSE, cluster_cols = FALSE,  angle_col = 0,
         legend = TRUE, legend_labels = c("None", "HR only", "RA only", "Both"),
         legend_breaks = c(0,1,2,3))

# network for all diet items
# matrix of pest occurrences
allSacc <- c(unique(insectData$`paste(insectData$order, insectData$sacc, sep = "")`))

dietMatrix <- matrix(ncol = 0, nrow = length(allSacc))
for (name in names(insectTaxonomy[36:39])){
  yesNo <- c()
  dietDF <- insectTaxonomy[[name]][-c(which(insectTaxonomy[[name]]$order=="NA")),]

  for (all in allSacc) {
    ifelse(all %in% unique(dietDF$`paste(insectData$order, insectData$sacc, sep = "")`), 
           yesNo <- append(yesNo, 1),
           yesNo <- append(yesNo, 0))
  }
dietMatrix <- cbind(dietMatrix, yesNo)
}

rownames(dietMatrix) <- c(allSacc)
colnames(dietMatrix) <- c(names(insectTaxonomy[36:39]))

dietMatrixCut <- filter_rows(dietMatrix)
dietMatrix <- dietMatrix[-c(dietMatrixCut),]

rownames(dietMatrix) <- c(abbreviate(c(rownames(dietMatrix)), minlength = 3, 
                                     strict = TRUE, method = "left.keep", use.classes = FALSE))

# diet networks
orderAbbr <- c(unique(rownames(dietMatrix)))
dietTotal <- list()
for (diet in orderAbbr) {
  dietDF <- data.frame(matrix(dietMatrix[which(str_detect(diet,rownames(dietMatrix))==TRUE),], ncol = 4))
  
  dietList <- c()
  for (col in 1:ncol(dietDF)) {
    dietList <- append(dietList, sum(dietDF[,col])) 
  }
  dietTotal <- append(dietTotal, list(dietList))
}

library(bipartite)


dietBipartite <- matrix(c(dietTotal[[1]], dietTotal[[2]], 
                          dietTotal[[3]], dietTotal[[4]],
                          dietTotal[[5]], dietTotal[[6]],
                          dietTotal[[7]], dietTotal[[8]]), ncol=8)
rownames(dietBipartite) <- c(LocationSeasonList[-c(5:6)])
colnames(dietBipartite) <- c(orderAbbr)

plotweb(dietBipartite, 
        text.rot=90, labsize = 1.5,  
        y.width.low=0.05,y.width.high=0.05,
        col.high="light blue",col.low="light green",
        col.interaction = "grey90", ybig = 1.5,
        x.lim = c(0,1.2), sequence = list(seq.low = c(LocationSeasonList[-c(5:6)]),
                                          seq.high = c(orderAbbr)), y.lim = c(-0.5,2)) 


visweb(dietBipartite, labsize = 0.3) 

group.venn(list(H.ruber = c(unique(HRorder$`paste(insectData$order, insectData$sacc, sep = "")`)),
                R.alcyone = c(unique(RAorder$`paste(insectData$order, insectData$sacc, sep = "")`))),
           label = FALSE, lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(220, 150), cat.dist = 0.05, fill = c("lightblue2", "lightpink"))


for (name in names(insectTaxonomy[30:41])) {
  insectTaxonomy <- append(insectTaxonomy, 
                           list(insectTaxonomy[[name]][-c(which(insectTaxonomy[[name]]$order=="NA")),]))
  
}

names(insectTaxonomy) <- c(names(insectTaxonomy[1:61]), "HR_KW_order", "HR_KD_order",
                                 "HR_AW_order", "HR_AD_order", "HR_BW_order",
                                 "HR_BD_order", "RA_KW_order", "RA_KD_order",
                                 "RA_AW_order", "RA_AD_order", "RA_BW_order",
                                 "RA_BD_order")

JessVenn(list(HRwet = unique(c(insectTaxonomy$HR_BW_order$`paste(insectData$order, insectData$sacc, sep = "")`)),
                               #insectTaxonomy$HR_AW_order$`paste(insectData$order, insectData$sacc, sep = "")`,
                               #insectTaxonomy$HR_BW_order$`paste(insectData$order, insectData$sacc, sep = "")`)),
              HRdry = unique(c(insectTaxonomy$HR_BD_order$`paste(insectData$order, insectData$sacc, sep = "")`)),
                               #insectTaxonomy$HR_AD_order$`paste(insectData$order, insectData$sacc, sep = "")`,
                               #insectTaxonomy$HR_BD_order$`paste(insectData$order, insectData$sacc, sep = "")`))),
              RAwet = unique(c(insectTaxonomy$RA_BW_order$`paste(insectData$order, insectData$sacc, sep = "")`)),
                               #insectTaxonomy$RA_AW_order$`paste(insectData$order, insectData$sacc, sep = "")`,
                               #insectTaxonomy$RA_BW_order$`paste(insectData$order, insectData$sacc, sep = "")`)),
              RAdry = unique(c(insectTaxonomy$RA_BD_order$`paste(insectData$order, insectData$sacc, sep = "")`))),
                               #insectTaxonomy$RA_AD_order$`paste(insectData$order, insectData$sacc, sep = "")`,
                               #insectTaxonomy$RA_BD_order$`paste(insectData$order, insectData$sacc, sep = "")`))),
         label = FALSE, lab.cex=1, cex = 1.5,
         cat.cex = 1.5, cat.dist = c(0.23,0.23, 0.12, 0.12))

      
# save focal datasets
writexl::write_xlsx(insectTaxonomy$HR_occurrences_OTUs, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\HR.data.xlsx")
writexl::write_xlsx(insectTaxonomy$RA_occurrences_OTUs, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\RA.data.xlsx")

# calculate Shannon's and Levin's index on wPOO and RRA


wPOO_Levins <- c()
for (name in names(insectTaxonomy[30:39])) {
  dietProp <- c()
  for (row in 1:nrow(insectTaxonomy[[name]])) {
    dietProp <- append(dietProp, 
                       (insectTaxonomy[[name]][row, "wPOOlist[1:nrow(insectTaxonomy[[name]])]"]/sum(insectTaxonomy[[name]]$`wPOOlist[1:nrow(insectTaxonomy[[name]])]`))^2)
  }
  Levins <- 1/sum(dietProp)
  wPOO_Levins <- append(wPOO_Levins, (1/(length(insectTaxonomy[[name]][,1])-1))*(Levins-1))
}


RRA_Levins <- c()
for (name in names(insectTaxonomy[30:39])) {
  dietProp <- c()
  for (row in 1:nrow(insectTaxonomy[[name]])) {
    dietProp <- append(dietProp, 
                       (insectTaxonomy[[name]][row, "RRAlist[c(1:nrow(insectTaxonomy[[name]]))]"]/sum(insectTaxonomy[[name]]$`RRAlist[c(1:nrow(insectTaxonomy[[name]]))]`))^2)
  }
  Levins <- 1/sum(dietProp)
  RRA_Levins <- append(RRA_Levins, (1/(length(insectTaxonomy[[name]][,1])-1))*(Levins-1))
}

HRlevinsData <- data.frame(LS = rep(LocationSeasonList, times = 2),
                         Metric = c(rep("wPOO", times = 6), rep("RRA", times = 6)),
                         Levins = c(wPOO_Levins[1:6], RRA_Levins[1:6])) 

HRlevinsData$LS <- factor(HRlevinsData$LS,
                        levels = c(LocationSeasonList))

y_expression <- expression("Levin's index" ~ "" ~(B[A]))

HRlevinsPlot <- ggplot(HRlevinsData, aes(x = LS, y = Levins, fill = Metric)) + 
  geom_bar(stat="identity", position = "dodge") +
  theme_minimal() + scale_fill_manual(values = c("darkolivegreen3", "darkorange")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  labs(x = "Location - Season", y = y_expression)

compareShannonDF <- compareShannonStats[compareShannonStats$LS=="HR - Bokito - wet",]
wilcox.test(compareShannonDF$Counts, mu = RRA_Shannon[5])

# plot difference in Shannon's and Levin's indices

library(ggplot2)
library(reshape2)
library(dplyr)

df <- insectTaxonomy$RA_AD_occurrence
for (row in 1:nrow(df)) {
  for (col in 1:ncol(df)) {
    ifelse(df[row,col]==1,
           df[row,col] <- insectTaxonomy$RA_AD_info[row,"RRAlist[c(1:nrow(insectTaxonomy[[name]]))]"],
           df[row,col] <- df[row,col])
  }
}


insectTaxonomy <- append(insectTaxonomy, list(HR_KW_RRA = HR_KW_RRAlevins,
                                              HR_KD_RRA = HR_KD_RRAlevins,
                                              HR_AW_RRA = HR_AW_RRAlevins,
                                              HR_AD_RRA = HR_AD_RRAlevins,
                                              HR_BW_RRA = HR_BW_RRAlevins,
                                              HR_BD_RRA = HR_BD_RRAlevins,
                                              RA_KW_RRA = RA_KW_RRAlevins,
                                              RA_KD_RRA = RA_KD_RRAlevins,
                                              RA_AW_RRA = RA_AW_RRAlevins, 
                                              RA_AD_RRA = RA_AD_RRAlevins))

allRRAshannon <- c()
for (name in names(insectTaxonomy[52:61])) {
  allRRAshannon <- append(allRRAshannon, ShannonIndexFunction(insectTaxonomy[[name]]))
}

RRAshannonStats <- data.frame(LS = c(rep("HR - Konye - wet", times = length(insectTaxonomy$HR_KW_RRA)),
                                        rep("HR - Konye - dry", times = length(insectTaxonomy$HR_KD_RRA)),
                                        rep("HR - Ayos - wet", times = length(insectTaxonomy$HR_AW_RRA)),
                                        rep("HR - Ayos - dry", times = length(insectTaxonomy$HR_AD_RRA)),
                                        rep("HR - Bokito - wet", times = length(insectTaxonomy$HR_BW_RRA)),
                                        rep("HR - Bokito - dry", times = length(insectTaxonomy$HR_BD_RRA)),
                                        rep("RA - Konye - wet", times = length(insectTaxonomy$RA_KW_RRA)),
                                        rep("RA - Konye - dry", times = length(insectTaxonomy$RA_KD_RRA)),
                                        rep("RA - Ayos - wet", times = length(insectTaxonomy$RA_AW_RRA)),
                                        rep("RA - Ayos - dry", times = length(insectTaxonomy$RA_AD_RRA))),
                                 Counts = c(allRRAshannon))


RRAshannonMeanList <- c()
for (LS in unique(RRAshannonStats$LS)) {
  Sdf <- RRAshannonStats[RRAshannonStats$LS==LS,]
  RRAshannonMeanList <- append(RRAshannonMeanList, mean(Sdf$Counts))
}

incidenceLevins <- data.frame(Species = c(rep("H. ruber", times = 12), rep("R. alcyone", times = 8)),
                              LS = c(rep(LocationSeasonList, times = 2), rep(LocationSeasonList[-c(5:6)], times = 2)),
                              Level = (c(rep("individual wPOO", times = 6), 
                                        rep("population wPOO", times = 6),
                                        rep("individual wPOO", times = 4),
                                        rep("population wPOO", times = 4))),
                              Levins = c(wPOOlevinsMeanList[1:6], wPOO_Levins[1:6],
                                          wPOOlevinsMeanList[7:10], wPOO_Levins[7:10]))

incidenceLevins$LS <- factor(incidenceLevins$LS,
                             levels = c("Konye - wet", "Konye - dry",
                                        "Ayos - wet", "Ayos - dry",
                                        "Bokito - wet", "Bokito - dry"))

g1 <- ggplot(incidenceLevins, aes(x=LS, y = Levins, color = Level, group = LS)) + 
  geom_point(aes(color = Level),size=4) +
  geom_path((aes(x=LS, y=Levins, group=LS)), arrow=arrow(type = "closed", 
                                                         length = unit(0.2, "cm")),
                                                         #ends = "both"), 
            color = "black", size = 0.5) +
  theme_minimal() + ylim(0,1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(x = "Location - Season", y = "Shannon Index") +
  scale_colour_manual(values = c("gold2", "slateblue")) +
  facet_wrap(vars(Species))

g1
