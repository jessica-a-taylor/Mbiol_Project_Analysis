# Function for commenting out code

#####
FormatComment <- function() {     
  y <- as.list(readClipboard())     
  spacer <- function(x) 
    paste("#", paste("   ", collapse=""), x, sep="")     
  z <- sapply(y, spacer)    
  zz <- as.matrix(as.data.frame(z))     
  dimnames(zz) <- list(c(rep("", nrow(zz))), c(""))     
  writeClipboard(noquote(zz), format = 1)     
  return(noquote(zz)) }
#####

# import seq data
seq.data <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_consensus (1).xlsx")

# filter for insect OTUs
seq.insecta <- which(seq.data$class=="Insecta")
insect.data <- seq.data[c(seq.insecta),]
insect.data <- insect.data[,-c(56,185,272,274:291)]

# removing OTUs with less than 10 reads
# function to remove rows where all cells are below a certain threshold
#####
#   filter.rows <- function(x, y, z) {                  
#     z <- c()                                          
#     for (row in 1:nrow(x)) {                          
#       which.row <- which(x[row,] < y)                 
#       if (length(which.row)==length(x[row,])) {       
#         z <- append(z, row)                           
#       }                                               
#     }                                                 
#     return(z)                                         
#   }
#####     
   # insect.data.cut <- filter.rows(insect.data[,2:270], 10, z)
                                                                   
   # insect.data <- insect.data[-c(insect.data.cut),]                
                                
# normalise to remove OTUs that are represented by less than 1% of reads in each sample
   # new.insect.data <- as.data.frame(insect.data[,-c(1,271:276)])
#####
#   sum.reads <- c()                                                   
#      for (nid.sum in 1:ncol(new.insect.data)) {                      
#        sum.reads <- append(sum.reads, sum(new.insect.data[,nid.sum]))
#      }                                                               
#                                                                      
#      new.insect.data <- new.insect.data[,-c(247,248)]                
#      sum.reads <- c()                                                
#      for (nid.sum in 1:ncol(new.insect.data)) {                      
#        sum.reads <- append(sum.reads, sum(new.insect.data[,nid.sum]))
#      }
#                                                                      
#      colnames(new.insect.data) <- c(sum.reads)                       
#                                                                      
#      for (idc in as.character(sum.reads)) {                          
#        normalised.diet <- (new.insect.data[,idc]/as.numeric(idc))*100
#        new.insect.data[,idc] <- normalised.diet                      
#      }
#####
   
#   new.insect.data.cut <- filter.rows(new.insect.data, 1, z)                                                                 
#   insect.data <- insect.data[-c(new.insect.data.cut),]                                                                      

insect.data <- insect.data[,-c(248,249)]

# diet calculations
new.insect.data <- as.data.frame(insect.data[,-c(1,269:275)])

# function for calculating FOO
#####
FOO.function <- function(df) {                                                                
  insect.FOO <- c()                                                                 
  for (seq_row in 1:nrow(df)) {                                                
    I <- 0                                                                                  
    for (seq_col in 1:ncol(df)) {                                              
      ifelse(df[seq_row, seq_col] > 0,                                         
             I <- I + 1, I <- I) 
    }                                                                                       
    insect.FOO <- append(insect.FOO, (1/length(names(df)))*I*100)
  }  
  return(insect.FOO)
}
#####

insect.FOO <- FOO.function(new.insect.data)
insect.data <- cbind(insect.data, insect.FOO)

# function for calculating POO
#####
POO.function <- function(df) {                                      
  insect.POO.counts <- c()                                        
  for (seq_row in 1:nrow(df)) {                      
    I <- 0                                                        
    for (seq_col in 1:ncol(df)) {                    
      ifelse(df[seq_row, seq_col] > 0,               
             I <- I + 1, I <- I)                                  
    }  
    insect.POO.counts <- append(insect.POO.counts, I)
  }                                                             
  sum.insect.POO.counts <- sum(insect.POO.counts)                 
  insect.POO <- c()                                               
  for (i in insect.POO.counts) {                                  
    insect.POO <- append(insect.POO,(i/sum.insect.POO.counts)*100)
  }    
  return(insect.POO)
}
#####

insect.POO <- POO.function(new.insect.data)
insect.data <- cbind(insect.data, insect.POO)

# function for calculating wPOO
#####
   wPOO.function <- function(df) {                                                                              
     for (row.no in 1:nrow(df)) {                                                                
       for (col.no in 1:ncol(df)) {                                                              
         if (df[row.no, col.no] > 0) {                                                            
           df[row.no, col.no] <- 1
         }                                                                                                    
       }                                                                                                      
     }
                                                                                                              
     insect.sample.diet <- c()                                                                                
     for (seq_col in 1:ncol(df)) {                                                               
     J <- 0
       for (seq_row in 1:nrow(df)) {                                                             
         ifelse(df[seq_row, seq_col] > 0, J <- J + 1, J <- J)                                    
       }                                                                                                      
       insect.sample.diet <- append(insect.sample.diet, J)                                                    
     }
     colnames(df) <- c(insect.sample.diet)                                                       

     df <- as.data.frame(df)                                                        
     reads <- c()
     for (insect.col in 1:ncol(df)) {                                                            
       reads <- append(reads, df[,insect.col])                                                   
     }                                                                                                        

     insect.df <- data.frame(Row = rep(1:length(df[,1]), times = length(df[1,])),                                                   
                             Reads = reads,                                                                   
                             Total.reads = rep(insect.sample.diet, each = length(df[,1])))                               

     insect.df <- cbind(insect.df, insect.df$Reads/insect.df$Total.reads)                                     
                                                                                                              
     insect.value.df <- df                                                                       
     for (row in 1:nrow(insect.value.df)) {
       new.insect.df <- insect.df[insect.df$Row==row,]                                                        
       insect.value.df[row,] <- new.insect.df$`insect.df$Reads/insect.df$Total.reads`                         
     }                                                                                                        

     insect.wPOO <- c()                                                                                       
     for (value_row in 1:nrow(insect.value.df)) {                                                             
       insect.wPOO <- append(insect.wPOO, 100*(1/length(insect.sample.diet))*sum(insect.value.df[value_row,]))
     }  
     return(insect.wPOO)
   }
#####

insect.wPOO <- wPOO.function(new.insect.data)
insect.data <- cbind(insect.data, insect.wPOO)

writexl::write_xlsx(insect.data, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\insect.filtered.data.xlsx")

# import habitat data
# select columns of interest and remove NAs
#####
farms <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\FecesSampleDatabase6Nov2019_MISEQ1&MISEQ2.xlsx", sheet = "Samples")                        
farms <- data.frame(farms$Lab.nbr...26, farms$animal, farms$Species, farms$Location,
                    farms$Site, as.character.POSIXt(farms$Date))                     
farms <- farms[-145,]                                                                
                                                                                     
farms.na <- which(is.na(farms$farms.Lab.nbr...26), arr.ind=TRUE)                     
farms <- farms[-c(farms.na),]

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
farms.ln <- as.character(farms$farms.Lab.nbr...26)
zbj_names <- names(insect.data[,2:268])           
                                                  
Lab.no <- c()                                     
for (n in zbj_names) {                            
for (l in farms.ln) {                           
  if (l==n)
    Lab.no <- append(Lab.no, as.numeric(n))     
  }                                               
}
#####

# compile information of interest from both datasets
#####
library(stringr)
zbj_df <- data.frame()                                                                                                                  
                                                                                                                                           
for (i in Lab.no) {                                                                                                                     
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
library(writexl)
write_xlsx(zbj_final, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_final.xlsx")

insect.data <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\insect.filtered.data.xlsx")
zbj_final <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_final.xlsx")

# create hash for matching lab numbers to species name
#####
predator.names <- hash()                             
for (row in 1:nrow(zbj_final)) {                     
  pred.nbr <- as.character(zbj_final[row, "Lab.nbr"])
  pred.sp <- as.character(zbj_final[row,"Species"])  
  predator.names[[pred.nbr]] <- pred.sp              
}
#####

# use hash to create a matrix of sequence data for the focal species
#####
zbj.insect.data <- insect.data[,c(2:268)]                                                  
zbj.insect.data <- zbj.insect.data[,-c(218, 222)]                                     
new.zbj.insect.data <- zbj.insect.data                                                     

pred.colnames <- c()                                                                       
for (name in names(new.zbj.insect.data)) {                                                 
  pred.colnames <- append(pred.colnames, predator.names[[name]])                           
}

colnames(new.zbj.insect.data) <- c(pred.colnames)                                          
                                                                                           
HR.focal.nbrs <- which(names(new.zbj.insect.data) == "Hipposideros ruber", arr.ind = TRUE) 
RA.focal.nbrs <- which(names(new.zbj.insect.data) == "Rhinolophus alcyone", arr.ind = TRUE)
                                                                                              
all.focal.nbrs <- c(HR.focal.nbrs, RA.focal.nbrs)                                          
all.focal.nbrs <- sort(all.focal.nbrs)                                                     
                                                                                           
zbj.insect.data.names <- names(zbj.insect.data)
focal.nbrs <- zbj.insect.data.names[all.focal.nbrs]                                        
zbj.focal <- zbj.insect.data[,c(focal.nbrs)]                                               
#####

# function for converting seq data to occurrence data
##### 
seq.to.occurrence <- function(df) { 
  ifelse(df > 0,                    
         df <- 1, df <- 0)                 
}
#####

# function for converting seq data to RRA data
seq.to.RRA <- function(df) {
  new.df <- df
  for (row in 1:nrow(df)) {
    for (col in 1:ncol(df)) {
      if (df[row,col] >= 1) {
        new.df[row,col] <- (df[row,col]/sum(df[,col]))*100
      }
    }
  }
  return(new.df)
}

zbj.focal.RRA <- seq.to.RRA(zbj.focal)

# flip columns and rows
zbj.focal <- t(zbj.focal)
zbj.focal <- as.data.frame(zbj.focal)

zbj.focal.RRA <- t(zbj.focal.RRA)
zbj.focal.RRA <- as.data.frame(zbj.focal.RRA)

library(dplyr)

zbj.focal <- data.frame(lapply(zbj.focal,seq.to.occurrence))

# add columns with sample and habitat info
#####                                                                                          
new.zbj.insect.data.names <- names(new.zbj.insect.data)                                    
focal.names <- new.zbj.insect.data.names[all.focal.nbrs]
zbj.focal <- as.data.frame(zbj.focal)                       
zbj.focal <- cbind(zbj.focal, focal.nbrs)                   
zbj.focal <- cbind(zbj.focal, focal.names)                  
                                                              
focal.info <- data.frame()                                  
for (nbr in as.numeric(focal.nbrs)) {                       
  for (row in 1:nrow(zbj_final)) {                          
    ifelse(zbj_final[row, "Lab.nbr"] == nbr,
           focal.info <- rbind(focal.info, zbj_final[row,]),
           focal.info <- focal.info)                        
  }                                                         
}                                                           
                                                              
zbj.focal <- cbind(zbj.focal, focal.info[,c(4,6)])
colnames(zbj.focal) <- c(insect.data$OTU, "Lab.nbr", "Species", "Location", "Season")

zbj.focal.RRA <- cbind(zbj.focal.RRA, focal.nbrs)                   
zbj.focal.RRA <- cbind(zbj.focal.RRA, focal.names)  
zbj.focal.RRA <- cbind(zbj.focal.RRA, focal.info[,c(4,6)])
colnames(zbj.focal.RRA) <- c(insect.data$OTU, "Lab.nbr", "Species", "Location", "Season")

focal.info <- focal.info[,-c(1,2,5)]

#####
library(writexl)
write_xlsx(zbj.focal, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\focal.species.data.xlsx")
write_xlsx(zbj.focal.RRA, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\focal.species.data.RRA.xlsx")

#####
filter.cols <- function(x, y, z) {                 
  z <- c()                                          
  for (col in 1:ncol(x)) {                         
    which.col <- which(x[,col] < y)                
    if (length(which.col)==length(x[,col])) {      
      z <- append(z, col) 
    }                                               
  } 
  return(z)
}

zbj.focal.RRA.cut <- filter.cols(zbj.focal.RRA[1:13979], 0.000001, z)
new.zbj.rra <- zbj.focal.RRA[,-c(zbj.focal.RRA.cut)]

zbj.focal.cut <- filter.cols(zbj.focal[1:13979], 1, z)
new.zbj.focal <- zbj.focal[,-c(zbj.focal.cut)]

# total number of reads for all species == OTU.sum
#####
insect.data.reads <- insect.data[,c(2:268)]          
OTU.all <- c()                                            
for (col in 1:ncol(insect.data.reads)) {                  
  OTU.all <- append(OTU.all, sum(insect.data.reads[,col]))
}                                                         
                                                             
OTU.sum <- sum(OTU.all)
#####

# isolate sequence data for focal species 
#####
HR.nbrs <- zbj.insect.data.names[HR.focal.nbrs]
HR.focal.df <- zbj.insect.data[,c(HR.nbrs)]
filter.HR.focal.df <- filter.rows(HR.focal.df[1:79], 1, z)
HR.focal.df <- HR.focal.df[-c(filter.HR.focal.df),]

RA.nbrs <- zbj.insect.data.names[RA.focal.nbrs]
RA.focal.df <- zbj.insect.data[,c(RA.nbrs)]
filter.RA.focal.df <- filter.rows(RA.focal.df[1:82], 1, z)
RA.focal.df <- RA.focal.df[-c(filter.RA.focal.df),]
#####

# number of reads for H.ruber
#####
HR.OTU.all <- c()
for (col in 1:ncol(HR.focal.df[1:79])) {
  HR.OTU.all <- append(HR.OTU.all, sum(HR.focal.df[,col]))
}

HR.OTU.sum <- sum(HR.OTU.all)
#####

# number of reads for R.alcyone
#####
RA.OTU.all <- c()
for (col in 1:ncol(RA.focal.df[1:82])) {
  RA.OTU.all <- append(RA.OTU.all, sum(RA.focal.df[,col]))
}

RA.OTU.sum <- sum(RA.OTU.all)
#####

# taxonomic resolution
#####
insect.data.order <- insect.data[,c(1,271:274,276)]

HR.focal.df <- cbind(HR.focal.df, insect.data.order[-(filter.HR.focal.df),])

RA.focal.df <- cbind(RA.focal.df, insect.data.order[-(filter.RA.focal.df),])


HRsp <- HR.focal.df$species
HR.no.sp <- HRsp[-c(which(str_detect(HRsp, "sp.")))]

is.na <- c()
for (row in 1:length(HR.no.sp)) {
  ifelse(HR.no.sp[row]=="NA",
         is.na <- append(is.na, row),
         is.na <- is.na)
}

(205/length(HR.focal.df$order))*100

new.HR.focal.df <- as.data.frame(HR.focal.df)
new.RA.focal.df <- as.data.frame(RA.focal.df)
#####


# proportion of Lepidoptera
#####
RA.order <- RA.focal.df[-c(which(RA.focal.df$order=="NA")),]
RA_Lepidoptera <- length(which(RA.order$order=="Lepidoptera"))
(RA_Lepidoptera/length(RA.order$order))*100

length(unique(RA.order$family))
#####

HR.focal.df <- data.frame(lapply(HR.focal.df[,1:79],seq.to.occurrence))
colnames(HR.focal.df) <- c(HR.nbrs)

HR.focal.df.cut <- filter.rows(HR.focal.df, 1, z)

RA.focal.df <- data.frame(lapply(RA.focal.df[,1:82],seq.to.occurrence))
colnames(RA.focal.df) <- c(RA.nbrs)

RA.focal.df.cut <- filter.rows(RA.focal.df, 1, z)

new.HR.focal.df <- seq.to.RRA(HR.focal.df[,1:79])
new.HR.focal.df <- cbind(new.HR.focal.df, HR.focal.df[80:85])
colnames(new.HR.focal.df) <- c(HR.nbrs, "OTU", "order", "family", "genus", "species", "insect.wPOO")

new.RA.focal.df <- seq.to.RRA(RA.focal.df[,1:82])
new.RA.focal.df <- cbind(new.RA.focal.df, RA.focal.df[83:88])
colnames(new.RA.focal.df) <- c(RA.nbrs, "OTU", "order", "family", "genus", "species", "insect.wPOO")


# presence/absence data for each location-season
#####
new.HR.zbj.data <- HR.zbj.data
new.HR.zbj.data <- cbind(new.HR.zbj.data, paste(new.HR.zbj.data$Location, "-" ,new.HR.zbj.data$Season))
colnames(new.HR.zbj.data) <- c(names(HR.zbj.data), "Location - Season")

HR.KW.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Konye - wet",]
HR.KW.nbrs <- HR.KW.data$Lab.nbr

HR.KW.focal.df <- HR.focal.df[,c(c(as.character(HR.KW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
HR.KW.focal.trim <- filter.rows(HR.KW.focal.df[,1:13], 1, z)
HR.KW.focal.cut <- HR.KW.focal.df[-c(HR.KW.focal.trim),]

HR.KD.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Konye - dry",]
HR.KD.nbrs <- HR.KD.data$Lab.nbr

HR.KD.focal.df <- HR.focal.df[,c(c(as.character(HR.KD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
HR.KD.focal.trim <- filter.rows(HR.KD.focal.df[,1:24], 1, z)
HR.KD.focal.cut <- HR.KD.focal.df[-c(HR.KD.focal.trim),]

HR.AW.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Ayos - wet",]
HR.AW.nbrs <- HR.AW.data$Lab.nbr

HR.AW.focal.df <- HR.focal.df[,c(c(as.character(HR.AW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
HR.AW.focal.trim <- filter.rows(HR.AW.focal.df[,1:10], 1, z)
HR.AW.focal.cut <- HR.AW.focal.df[-c(HR.AW.focal.trim),]

HR.AD.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Ayos - dry",]
HR.AD.nbrs <- HR.AD.data$Lab.nbr

HR.AD.focal.df <- HR.focal.df[,c(c(as.character(HR.AD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
HR.AD.focal.trim <- filter.rows(HR.AD.focal.df[,1:25], 1, z)
HR.AD.focal.cut <- HR.AD.focal.df[-c(HR.AD.focal.trim),]

HR.BW.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Bokito - wet",]
HR.BW.nbrs <- HR.BW.data$Lab.nbr

HR.BW.focal.df <- HR.focal.df[,c(c(as.character(HR.BW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
HR.BW.focal.trim <- filter.rows(HR.BW.focal.df[,1:4], 1, z)
HR.BW.focal.cut <- HR.BW.focal.df[-c(HR.BW.focal.trim),]

HR.BD.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Bokito - dry",]
HR.BD.nbrs <- HR.BD.data$Lab.nbr

HR.BD.focal.df <- HR.focal.df[,c(c(as.character(HR.BD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
HR.BD.focal.trim <- filter.rows(HR.BD.focal.df[,1:3], 1, z)
HR.BD.focal.cut <- HR.BD.focal.df[-c(HR.BD.focal.trim),]

new.RA.zbj.data <- RA.zbj.data
new.RA.zbj.data <- cbind(new.RA.zbj.data, paste(new.RA.zbj.data$Location, "-" ,new.RA.zbj.data$Season))
colnames(new.RA.zbj.data) <- c(names(RA.zbj.data), "Location - Season")

RA.KW.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Konye - wet",]
RA.KW.nbrs <- RA.KW.data$Lab.nbr

RA.KW.focal.df <- RA.focal.df[,c(c(as.character(RA.KW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
RA.KW.focal.trim <- filter.rows(RA.KW.focal.df[,1:29], 1, z)
RA.KW.focal.cut <- RA.KW.focal.df[-c(RA.KW.focal.trim),]

RA.KD.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Konye - dry",]
RA.KD.nbrs <- RA.KD.data$Lab.nbr

RA.KD.focal.df <- RA.focal.df[,c(c(as.character(RA.KD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
RA.KD.focal.trim <- filter.rows(RA.KD.focal.df[,1:18], 1, z)
RA.KD.focal.cut <- RA.KD.focal.df[-c(RA.KD.focal.trim),]

RA.AW.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Ayos - wet",]
RA.AW.nbrs <- RA.AW.data$Lab.nbr

RA.AW.focal.df <- RA.focal.df[,c(c(as.character(RA.AW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
RA.AW.focal.trim <- filter.rows(RA.AW.focal.df[,1:20], 1, z)
RA.AW.focal.cut <- RA.AW.focal.df[-c(RA.AW.focal.trim),]

RA.AD.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Ayos - dry",]
RA.AD.nbrs <- RA.AD.data$Lab.nbr

RA.AD.focal.df <- RA.focal.df[,c(c(as.character(RA.AD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
RA.AD.focal.trim <- filter.rows(RA.AD.focal.df[,1:13], 1, z)
RA.AD.focal.cut <- RA.AD.focal.df[-c(RA.AD.focal.trim),]

RA.BW.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Bokito - wet",]
RA.BW.nbrs <- RA.BW.data$Lab.nbr

RA.BW.focal.df <- RA.focal.df[,c(c(as.character(RA.BW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
RA.BW.focal.trim <- which(RA.BW.focal.df==0)
RA.BW.focal.cut <- as.data.frame(RA.BW.focal.df[-c(RA.BW.focal.trim)])

RA.BD.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Bokito - dry",]
RA.BD.nbrs <- RA.BD.data$Lab.nbr

RA.BD.focal.df <- RA.focal.df[,c(c(as.character(RA.BD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
RA.BD.focal.trim <- which(RA.BD.focal.df==0)
RA.BD.focal.cut <- as.data.frame(RA.BD.focal.df[-c(RA.BD.focal.trim)])

#####

# RRA data for each location-season
#####
new.HR.KW.focal.df <- new.HR.focal.df[,c(c(as.character(HR.KW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.HR.KW.focal.trim <- filter.rows(new.HR.KW.focal.df[,1:13], 0.000001, z)
new.HR.KW.focal.cut <- new.HR.KW.focal.df[-c(new.HR.KW.focal.trim),]

new.HR.KD.focal.df <- new.HR.focal.df[,c(c(as.character(HR.KD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.HR.KD.focal.trim <- filter.rows(new.HR.KD.focal.df[,1:24], 0.000001, z)
new.HR.KD.focal.cut <- new.HR.KD.focal.df[-c(new.HR.KD.focal.trim),]

new.HR.AW.focal.df <- new.HR.focal.df[,c(c(as.character(HR.AW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.HR.AW.focal.trim <- filter.rows(new.HR.AW.focal.df[,1:10], 0.000001, z)
new.HR.AW.focal.cut <- new.HR.AW.focal.df[-c(new.HR.AW.focal.trim),]

new.HR.AD.focal.df <- new.HR.focal.df[,c(c(as.character(HR.AD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.HR.AD.focal.trim <- filter.rows(new.HR.AD.focal.df[,1:25], 0.000001, z)
new.HR.AD.focal.cut <- new.HR.AD.focal.df[-c(new.HR.AD.focal.trim),]

new.HR.BW.focal.df <- new.HR.focal.df[,c(c(as.character(HR.BW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.HR.BW.focal.trim <- filter.rows(new.HR.BW.focal.df[,1:4], 1, z)
new.HR.BW.focal.cut <- new.HR.BW.focal.df[-c(new.HR.BW.focal.trim),]

new.HR.BD.focal.df <- new.HR.focal.df[,c(c(as.character(HR.BD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.HR.BD.focal.trim <- filter.rows(new.HR.BD.focal.df[,1:3], 0.000001, z)
new.HR.BD.focal.cut <- new.HR.BD.focal.df[-c(new.HR.BD.focal.trim),]

new.RA.KW.focal.df <- new.RA.focal.df[,c(c(as.character(RA.KW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.RA.KW.focal.trim <- filter.rows(new.RA.KW.focal.df[,1:29], 1, z)
new.RA.KW.focal.cut <- new.RA.KW.focal.df[-c(new.RA.KW.focal.trim),]

new.RA.KD.focal.df <- new.RA.focal.df[,c(c(as.character(RA.KD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.RA.KD.focal.trim <- filter.rows(new.RA.KD.focal.df[,1:18], 0.000001, z)
new.RA.KD.focal.cut <- new.RA.KD.focal.df[-c(new.RA.KD.focal.trim),]

new.RA.AW.focal.df <- new.RA.focal.df[,c(c(as.character(RA.AW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.RA.AW.focal.trim <- filter.rows(new.RA.AW.focal.df[,1:20], 0.000001, z)
new.RA.AW.focal.cut <- new.RA.AW.focal.df[-c(new.RA.AW.focal.trim),]

new.RA.AD.focal.df <- new.RA.focal.df[,c(c(as.character(RA.AD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
new.RA.AD.focal.trim <- filter.rows(new.RA.AD.focal.df[,1:13], 0.000001, z)
new.RA.AD.focal.cut <- new.RA.AD.focal.df[-c(new.RA.AD.focal.trim),]

new.RA.BW.focal.df <- new.RA.focal.df[,c(c(as.character(RA.BW.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]

new.RA.BD.focal.df <- new.RA.focal.df[,c(c(as.character(RA.BD.nbrs)), "OTU", "order", "family", "genus", "species", "insect.wPOO")]
#####

# function to count OTUs per sample
#####
count.OTUs <- function(df) {                             
  HR.count.OTUs <- c()                                   
  for (col in 1:ncol(df)) {                              
    count <- which(df[,col]==1)                          
    HR.count.OTUs <- append(HR.count.OTUs, length(count))
  }                                                      
  return(HR.count.OTUs)                                  
  
} 
#####

# compare average number of OTUs between focal species
#####
HR.count.OTUs <- count.OTUs(HR.focal.df[,1:79])
RA.count.OTUs <- count.OTUs(RA.focal.df[,1:82])
library(ggplot2)
focal.compare.OTUs.df <- data.frame(Species = c(rep("H.ruber", times = length(HR.count.OTUs)),
                                                rep("R.alcyone", times = length(RA.count.OTUs))),
                                    Counts = c(HR.count.OTUs, RA.count.OTUs))

focal.compare.OTUs.boxplot <- ggplot(focal.compare.OTUs.df, aes(x = Species, y = Counts, fill = Species)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.text.y = element_text(size = 10)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Species", y = "Number of OTUs")

wilcox.test(HR.count.OTUs, RA.count.OTUs)
shapiro.test(focal.compare.OTUs.df$Counts)
#####

# number of OTUs in each Location-Season
HR.KW.OTUs <- count.OTUs(HR.KW.focal.cut[1:13])
HR.KD.OTUs <- count.OTUs(HR.KD.focal.cut[1:24])
HR.AW.OTUs <- count.OTUs(HR.AW.focal.cut[1:10])
HR.AD.OTUs <- count.OTUs(HR.AD.focal.cut[1:25])
HR.BW.OTUs <- count.OTUs(HR.BW.focal.cut[1:4])
HR.BD.OTUs <- count.OTUs(HR.BD.focal.cut[1:3])

RA.KW.OTUs <- count.OTUs(RA.KW.focal.cut[1:29])
RA.KD.OTUs <- count.OTUs(RA.KD.focal.cut[1:18])
RA.AW.OTUs <- count.OTUs(RA.AW.focal.cut[1:20])
RA.AD.OTUs <- count.OTUs(RA.AD.focal.cut[1:13])
RA.BW.OTUs <- length(which(RA.BW.focal.df[1]==1))
RA.BD.OTUs <- length(which(RA.BD.focal.df[1]==1))

focal.compare.OTUs.df <- data.frame(Species = rep("H.ruber", times = length(c(HR.KW.OTUs, HR.KD.OTUs,
                                                                                HR.AW.OTUs,HR.AD.OTUs,
                                                                                HR.BW.OTUs,HR.BD.OTUs))),
                                    LS = c(rep("Konye - wet", times = length(HR.KW.OTUs)),
                                           rep("Konye - dry", times = length(HR.KD.OTUs)),
                                           rep("Ayos - wet", times = length(HR.AW.OTUs)),
                                           rep("Ayos - dry", times = length(HR.AD.OTUs)),
                                           rep("Bokito - wet", times = length(HR.BW.OTUs)),
                                           rep("Bokito - dry", times = length(HR.BD.OTUs))),
                                    Counts = c(HR.KW.OTUs, HR.KD.OTUs,HR.AW.OTUs,HR.AD.OTUs,HR.BW.OTUs,HR.BD.OTUs))

focal.compare.OTUs.df$LS <- factor(focal.compare.OTUs.df$LS,
                                   levels = c("Konye - wet", "Konye - dry",
                                              "Ayos - wet", "Ayos - dry",
                                              "Bokito - wet", "Bokito - dry"))

focal.compare.OTUs.boxplot <- ggplot(focal.compare.OTUs.df, aes(x = LS, y = Counts, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 10),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
                                                axis.text.y = element_text(size = 10)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = "Number of OTUs") + ylim(20,55)

focal.compare.OTUs.boxplot + annotate("text",
                                      x = 1:length(table(focal.compare.OTUs.df$LS)),
                                      y = rep(53, times = 6),
                                      label = table(focal.compare.OTUs.df$LS),
                                      col = "red",
                                      vjust = - 1)


focal.compare.OTUs.df <- data.frame(LS = c(rep("HR - Konye - wet", times = length(HR.KW.OTUs)),
                                           rep("HR - Konye - dry", times = length(HR.KD.OTUs)),
                                           rep("HR - Ayos - wet", times = length(HR.AW.OTUs)),
                                           rep("HR - Ayos - dry", times = length(HR.AD.OTUs)),
                                           rep("HR - Bokito - wet", times = length(HR.BW.OTUs)),
                                           rep("HR - Bokito - dry", times = length(HR.BD.OTUs)),
                                           rep("RA - Konye - wet", times = length(RA.KW.OTUs)),
                                           rep("RA - Konye - dry", times = length(RA.KD.OTUs)),
                                           rep("RA - Ayos - wet", times = length(RA.AW.OTUs)),
                                           rep("RA - Ayos - dry", times = length(RA.AD.OTUs)),
                                           rep("RA - Bokito - wet", times = length(RA.BW.OTUs)),
                                           rep("RA - Bokito - dry", times = length(RA.BD.OTUs))),
                                    Counts = c(HR.KW.OTUs, HR.KD.OTUs,HR.AW.OTUs,HR.AD.OTUs,HR.BW.OTUs,HR.BD.OTUs,
                                               RA.KW.OTUs, RA.KD.OTUs,RA.AW.OTUs,RA.AD.OTUs,RA.BW.OTUs,RA.BD.OTUs))

compare.OTUs.kruskal <- kruskal.test(data = focal.compare.OTUs.df, x = focal.compare.OTUs.df$Counts,
                                     g = focal.compare.OTUs.df$LS, formula = Counts ~ LS)


# occurrence calculations for focal species
#####
HR.KW.wPOO <- wPOO.function(HR.KW.focal.cut[1:13])
HR.KD.wPOO <- wPOO.function(HR.KD.focal.cut[1:24])
HR.AW.wPOO <- wPOO.function(HR.AW.focal.cut[1:10])
HR.AD.wPOO <- wPOO.function(HR.AD.focal.cut[1:25])
HR.BW.wPOO <- wPOO.function(HR.BW.focal.cut[1:4])
HR.BD.wPOO <- wPOO.function(HR.BD.focal.cut[1:3])

RA.KW.wPOO <- wPOO.function(RA.KW.focal.cut[1:29])
RA.KD.wPOO <- wPOO.function(RA.KD.focal.cut[1:18])
RA.AW.wPOO <- wPOO.function(RA.AW.focal.cut[1:20])
RA.AD.wPOO <- wPOO.function(RA.AD.focal.cut[1:13])
RA.BW.wPOO <- "NA"
RA.BD.wPOO <- "NA"


# relative read abundance
# function for calculating RRA
#####
RRA.function <- function(df) {
  sum.read.counts <- c()
  
  for (row in 1:nrow(df)) {
    total.read.counts <- c()
    
    for (col in 1:ncol(df)) {
      
      total.read.counts <- append(total.read.counts, as.numeric(df[row,col]/sum(df[,col])))
    }
    sum.read.counts <- append(sum.read.counts, sum(total.read.counts))
  }
  return(sum.read.counts)
}
#####

HR.sum.read.counts <- RRA.function(new.HR.KW.focal.cut[1:13])
HR.KW.RRA <- (HR.sum.read.counts/length(ncol(new.HR.KW.focal.cut)))*100
HR.sum.read.counts <- RRA.function(new.HR.KD.focal.cut[1:24])
HR.KD.RRA <- (HR.sum.read.counts/length(ncol(new.HR.KD.focal.cut)))*100
HR.sum.read.counts <- RRA.function(new.HR.AW.focal.cut[1:10])
HR.AW.RRA <- (HR.sum.read.counts/length(ncol(new.HR.AW.focal.cut)))*100
HR.sum.read.counts <- RRA.function(new.HR.AD.focal.cut[1:25])
HR.AD.RRA <- (HR.sum.read.counts/length(ncol(new.HR.AD.focal.cut)))*100
HR.sum.read.counts <- RRA.function(new.HR.BW.focal.cut[1:4])
HR.BW.RRA <- (HR.sum.read.counts/length(ncol(new.HR.BW.focal.cut)))*100
HR.sum.read.counts <- RRA.function(new.HR.BD.focal.cut[1:3])
HR.BD.RRA <- (HR.sum.read.counts/length(ncol(new.HR.BD.focal.cut)))*100


RA.sum.read.counts <- RRA.function(new.RA.KW.focal.cut[1:29])
RA.KW.RRA <- (RA.sum.read.counts/length(ncol(new.RA.KW.focal.cut)))*100
RA.sum.read.counts <- RRA.function(new.RA.KD.focal.cut[1:18])
RA.KD.RRA <- (RA.sum.read.counts/length(ncol(new.RA.KD.focal.cut)))*100
RA.sum.read.counts <- RRA.function(new.RA.AW.focal.cut[1:20])
RA.AW.RRA <- (RA.sum.read.counts/length(ncol(new.RA.AW.focal.cut)))*100
RA.sum.read.counts <- RRA.function(new.RA.AD.focal.cut[1:13])
RA.AD.RRA <- (RA.sum.read.counts/length(ncol(new.RA.AD.focal.cut)))*100
RA.BW.RRA <- "NA"
RA.BD.RRA <- "NA"

insect.data.order <- insect.data[,c(1,271:274)]

HR.KW.focal.cut <- cbind(HR.KW.focal.cut, data.frame(HR.KW.wPOO, HR.KW.RRA))
colnames(HR.KW.focal.cut) <- c(HR.KW.nbrs, "OTU", "order", "family", "genus", "species",
                              "insect.wPOO","wPOO", "RRA")

HR.KD.focal.cut <- cbind(HR.KD.focal.cut, data.frame(HR.KD.wPOO, HR.KD.RRA))
colnames(HR.KD.focal.cut) <- c(HR.KD.nbrs, "OTU", "order", "family", "genus", "species",
                               "insect.wPOO","wPOO", "RRA")


HR.AW.focal.cut <- cbind(HR.AW.focal.cut, data.frame(HR.AW.wPOO, HR.AW.RRA))
colnames(HR.AW.focal.cut) <- c(HR.AW.nbrs, "OTU", "order", "family", "genus", "species",
                               "insect.wPOO","wPOO", "RRA")

HR.AD.focal.cut <- cbind(HR.AD.focal.cut, data.frame(HR.AD.wPOO, HR.AD.RRA))
colnames(HR.AD.focal.cut) <- c(HR.AD.nbrs, "OTU", "order", "family", "genus", "species",
                               "insect.wPOO","wPOO", "RRA")

HR.BW.focal.cut <- cbind(HR.BW.focal.cut, data.frame(HR.BW.wPOO, HR.BW.RRA))
colnames(HR.BW.focal.cut) <- c(HR.BW.nbrs, "OTU", "order", "family", "genus", "species",
                               "insect.wPOO","wPOO", "RRA")

HR.BD.focal.df <- cbind(HR.BD.focal.df, insect.data.order)
HR.BD.focal.df <- cbind(HR.BD.focal.df, data.frame(HR.BD.wPOO, HR.BD.RRA))
colnames(HR.BD.focal.df) <- c(HR.BD.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

RA.KW.focal.df <- cbind(RA.KW.focal.df, insect.data.order)
RA.KW.focal.df <- cbind(RA.KW.focal.df, data.frame(RA.KW.wPOO, RA.KW.RRA))
colnames(RA.KW.focal.df) <- c(RA.KW.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

RA.KD.focal.df <- cbind(RA.KD.focal.df, insect.data.order)
RA.KD.focal.df <- cbind(RA.KD.focal.df, data.frame(RA.KD.wPOO, RA.KD.RRA))
colnames(RA.KD.focal.df) <- c(RA.KD.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

RA.AW.focal.df <- cbind(RA.AW.focal.df, insect.data.order)
RA.AW.focal.df <- cbind(RA.AW.focal.df, data.frame(RA.AW.wPOO, RA.AW.RRA))
colnames(RA.AW.focal.df) <- c(RA.AW.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

RA.AD.focal.df <- cbind(RA.AD.focal.df, insect.data.order)
RA.AD.focal.df <- cbind(RA.AD.focal.df, data.frame(RA.AD.wPOO, RA.AD.RRA))
colnames(RA.AD.focal.df) <- c(RA.AD.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

RA.BW.focal.df <- cbind(RA.BW.focal.df, insect.data.order)
RA.BW.focal.df <- cbind(RA.BW.focal.df, data.frame(RA.BW.wPOO, RA.BW.RRA))
colnames(RA.BW.focal.df) <- c(RA.BW.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

RA.BD.focal.df <- cbind(RA.BD.focal.df, insect.data.order)
RA.BD.focal.df <- cbind(RA.BD.focal.df, data.frame(RA.BD.wPOO, RA.BD.RRA))
colnames(RA.BD.focal.df) <- c(RA.BD.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

new.HR.KW.focal.df <- cbind(new.HR.KW.focal.df, insect.data.order)
new.HR.KW.focal.df <- cbind(new.HR.KW.focal.df, data.frame(HR.KW.wPOO, HR.KW.RRA))
colnames(new.HR.KW.focal.df) <- c(HR.KW.nbrs, "OTU", "order", "family", "genus", "species",
                              "wPOO", "RRA")

new.HR.KD.focal.df <- cbind(new.HR.KD.focal.df, insect.data.order)
new.HR.KD.focal.df <- cbind(new.HR.KD.focal.df, data.frame(HR.KD.wPOO, HR.KD.RRA))
colnames(new.HR.KD.focal.df) <- c(HR.KD.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

new.HR.AW.focal.df <- cbind(new.HR.AW.focal.df, insect.data.order)
new.HR.AW.focal.df <- cbind(new.HR.AW.focal.df, data.frame(HR.AW.wPOO, HR.AW.RRA))
colnames(new.HR.AW.focal.df) <- c(HR.AW.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

new.HR.AD.focal.df <- cbind(new.HR.AD.focal.df, insect.data.order)
new.HR.AD.focal.df <- cbind(new.HR.AD.focal.df, data.frame(HR.AD.wPOO, HR.AD.RRA))
colnames(new.HR.AD.focal.df) <- c(HR.AD.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

new.RA.KW.focal.df <- cbind(new.RA.KW.focal.df, insect.data.order)
new.RA.KW.focal.df <- cbind(new.RA.KW.focal.df, data.frame(RA.KW.wPOO, RA.KW.RRA))
colnames(new.RA.KW.focal.df) <- c(RA.KW.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

new.RA.KD.focal.df <- cbind(new.RA.KD.focal.df, insect.data.order)
new.RA.KD.focal.df <- cbind(new.RA.KD.focal.df, data.frame(RA.KD.wPOO, RA.KD.RRA))
colnames(new.RA.KD.focal.df) <- c(RA.KD.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

new.RA.AW.focal.df <- cbind(new.RA.AW.focal.df, insect.data.order)
new.RA.AW.focal.df <- cbind(new.RA.AW.focal.df, data.frame(RA.AW.wPOO, RA.AW.RRA))
colnames(new.RA.AW.focal.df) <- c(RA.AW.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

new.RA.AD.focal.df <- cbind(new.RA.AD.focal.df, insect.data.order)
new.RA.AD.focal.df <- cbind(new.RA.AD.focal.df, data.frame(RA.AD.wPOO, RA.AD.RRA))
colnames(new.RA.AD.focal.df) <- c(RA.AD.nbrs, "OTU", "order", "family", "genus", "species",
                                  "wPOO", "RRA")

# diet composition - most common order
#####
not.na <- c()
for (row in 1:nrow(HR.focal.df)) {
  ifelse(HR.focal.df[row, "order"]!="NA",
         not.na <- append(not.na, row),
         not.na <- not.na)
}

HR.focal.order <- HR.focal.df[c(not.na),]
HR.max.wPOO <- max(HR.focal.order$insect.wPOO)
HR.top.order <- HR.focal.order[HR.focal.order$insect.wPOO==HR.max.wPOO,]
HR.top.order <- HR.top.order$order

HR.top.families <- HR.focal.order[HR.focal.order$order==HR.top.order,]
HR.top.families <- unique(HR.top.families$family)

# repeat for RRA data
not.na <- c()
for (row in 1:nrow(new.HR.focal.df)) {
  ifelse(new.HR.focal.df[row, "order"]!="NA",
         not.na <- append(not.na, row),
         not.na <- not.na)
}

new.HR.focal.order <- new.HR.focal.df[c(not.na),]
HR.max.RRA <- max(new.HR.focal.order$HR.focal.RRA)
new.HR.top.order <- new.HR.focal.order[new.HR.focal.order$HR.focal.RRA==HR.max.RRA,]
new.HR.top.order <- new.HR.top.order$order

new.HR.top.families <- new.HR.focal.order[new.HR.focal.order$order==new.HR.top.order,]
new.HR.top.families <- unique(new.HR.top.families$family)

#####

# Shannon diversity 
# function for calculating Shannon index
#####
library(vegan)
Shannon.index.function <- function(df) {                                 
  Shannon.index <- c()                                                   
  
  for (col in 1:ncol(df)) {                                               
    Shannon.index <- append(Shannon.index, diversity(df[,col], "shannon"))
  }                                                                       
  return(Shannon.index)                                                  
} 
#####

# Shannon for each location-season
#####
HR.KW.Shannon.index <- Shannon.index.function(HR.KW.focal.df[,1:13])
HR.KD.Shannon.index <- Shannon.index.function(HR.KD.focal.df[,1:24])
HR.AW.Shannon.index <- Shannon.index.function(HR.AW.focal.df[,1:10])
HR.AD.Shannon.index <- Shannon.index.function(HR.AD.focal.df[,1:25])
HR.BW.Shannon.index <- Shannon.index.function(HR.BW.focal.df[,1:4])
HR.BD.Shannon.index <- Shannon.index.function(HR.BD.focal.df[,1:3])

RA.KW.Shannon.index <- Shannon.index.function(RA.KW.focal.df[,1:29])
RA.KD.Shannon.index <- Shannon.index.function(RA.KD.focal.df[,1:18])
RA.AW.Shannon.index <- Shannon.index.function(RA.AW.focal.df[,1:20])
RA.AD.Shannon.index <- Shannon.index.function(RA.AD.focal.df[,1:13])
RA.BW.Shannon.index <- diversity(RA.BW.focal.df[,1], "shannon")
RA.BD.Shannon.index <- diversity(RA.BD.focal.df[,1], "shannon")

HR.compare.Shannon <- data.frame(Species = rep("H.ruber", times = length(c(RA.KW.Shannon.index, RA.KD.Shannon.index,
                                                                           RA.AW.Shannon.index, RA.AD.Shannon.index,
                                                                           RA.BW.Shannon.index, RA.BD.Shannon.index))),
                                 LS = c(rep("Konye - wet", times = length(RA.KW.Shannon.index)),
                                        rep("Konye - dry", times = length(RA.KD.Shannon.index)),
                                        rep("Ayos - wet", times = length(RA.AW.Shannon.index)),
                                        rep("Ayos - dry", times = length(RA.AD.Shannon.index)),
                                        rep("Bokito - wet", times = length(RA.BW.Shannon.index)),
                                        rep("Bokito - dry", times = length(RA.BD.Shannon.index))),
                                 Shannon.index = c(RA.KW.Shannon.index, RA.KD.Shannon.index, 
                                                   RA.AW.Shannon.index, RA.AD.Shannon.index,
                                                   RA.BW.Shannon.index, RA.BD.Shannon.index))

HR.compare.Shannon$LS <- factor(HR.compare.Shannon$LS,
                                   levels = c("Konye - wet", "Konye - dry",
                                              "Ayos - wet", "Ayos - dry",
                                              "Bokito - wet", "Bokito - dry"))

HR.compare.Shannon.boxplot <- ggplot(HR.compare.Shannon, aes(x = LS, y = Shannon.index, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = "Shannon index") + ylim(2.5,7.1)

HR.compare.Shannon.boxplot + annotate("text",
                                      x = 1:length(table(HR.compare.Shannon$LS)),
                                      y = rep(6.75, times = 6),
                                      label = table(HR.compare.Shannon$LS),
                                      col = "red",
                                      vjust = - 1)


focal.compare.Shannon.df <- data.frame(LS = c(rep("HR - Konye - wet", times = length(HR.KW.Shannon.index)),
                                           rep("HR - Konye - dry", times = length(HR.KD.Shannon.index)),
                                           rep("HR - Ayos - wet", times = length(HR.AW.Shannon.index)),
                                           rep("HR - Ayos - dry", times = length(HR.AD.Shannon.index)),
                                           rep("HR - Bokito - wet", times = length(HR.BW.Shannon.index)),
                                           rep("HR - Bokito - dry", times = length(HR.BD.Shannon.index)),
                                           rep("RA - Konye - wet", times = length(RA.KW.Shannon.index)),
                                           rep("RA - Konye - dry", times = length(RA.KD.Shannon.index)),
                                           rep("RA - Ayos - wet", times = length(RA.AW.Shannon.index)),
                                           rep("RA - Ayos - dry", times = length(RA.AD.Shannon.index)),
                                           rep("RA - Bokito - wet", times = length(RA.BW.Shannon.index)),
                                           rep("RA - Bokito - dry", times = length(RA.BD.Shannon.index))),
                                    Counts = c(new.HR.KW.Shannon.index, new.HR.KD.Shannon.index, new.HR.AW.Shannon.index,
                                               new.HR.AD.Shannon.index, new.HR.BW.Shannon.index, new.HR.BD.Shannon.index,
                                               new.RA.KW.Shannon.index, new.RA.KD.Shannon.index, new.RA.AW.Shannon.index,
                                               new.RA.AD.Shannon.index, new.RA.BW.Shannon.index, new.RA.BD.Shannon.index))

shapiro.test(focal.compare.Shannon.df$Counts)
compare.Shannon.kruskal <- kruskal.test(data = focal.compare.Shannon.df, x = focal.compare.Shannon.df$Counts,
                                     g = focal.compare.Shannon.df$LS, formula = Counts ~ LS)

#####

# Niche breadth
# function for calculating Levin's index
#####
Levins.function <- function(df) {
  sum.diet.prop <- c()
  for (col in 1:ncol(df)) {
    diet.prop <- c()
    for (row in 1:nrow(df)) {
      diet.prop <- append(diet.prop, (df[row,col]/sum(df[,col]))^2)
    }
    sum.diet.prop <- append(sum.diet.prop, sum(diet.prop))
  }
  
  Levins <- c()
  Levins.std <- c()
  for (prop in sum.diet.prop) {
    std <- 1/prop
    Levins <- append(Levins, std)
    Levins.std <- append(Levins.std, ((std-1)/(length(df[,1])-1)))
  }
  
  return(list(All.Levins = Levins,
              All.Levins.std = Levins.std,
              Mean.Levins = mean(Levins.std), 
              SD.Levins = sd(Levins.std)))
}
#####

# Levins for each location-season
#####
HR.KW.Levins.index <- Levins.function(HR.KW.focal.cut[,1:13])
HR.KD.Levins.index <- Levins.function(HR.KD.focal.cut[,1:24])
HR.AW.Levins.index <- Levins.function(HR.AW.focal.cut[,1:10])
HR.AD.Levins.index <- Levins.function(HR.AD.focal.cut[,1:25])
HR.BW.Levins.index <- Levins.function(HR.BW.focal.cut[,1:4])
HR.BD.Levins.index <- Levins.function(HR.BD.focal.cut[,1:3])

RA.KW.Levins.index <- Levins.function(RA.KW.focal.cut[,1:29])
RA.KD.Levins.index <- Levins.function(RA.KD.focal.cut[,1:18])
RA.AW.Levins.index <- Levins.function(RA.AW.focal.cut[,1:20])
RA.AD.Levins.index <- Levins.function(RA.AD.focal.cut[,1:13])

focal.compare.Levins.df <- data.frame(Species = rep("H.ruber", times = length(c(HR.KW.Levins.index$All.Levins.std, HR.KD.Levins.index$All.Levins.std,
                                                                          HR.AW.Levins.index$All.Levins.std, HR.AD.Levins.index$All.Levins.std,
                                                                          HR.BW.Levins.index$All.Levins.std, HR.BD.Levins.index$All.Levins.std))),
                                LS = c(rep("Konye - wet", times = length(HR.KW.Levins.index$All.Levins.std)),
                                       rep("Konye - dry", times = length(HR.KD.Levins.index$All.Levins.std)),
                                       rep("Ayos - wet", times = length(HR.AW.Levins.index$All.Levins.std)),
                                       rep("Ayos - dry", times = length(HR.AD.Levins.index$All.Levins.std)),
                                       rep("Bokito - wet", times = length(HR.BW.Levins.index$All.Levins.std)),
                                       rep("Bokito - dry", times = length(HR.BD.Levins.index$All.Levins.std))),
                                Levins.index = c(HR.KW.Levins.index$All.Levins.std, HR.KD.Levins.index$All.Levins.std, 
                                                 HR.AW.Levins.index$All.Levins.std, HR.AD.Levins.index$All.Levins.std,
                                                 HR.BW.Levins.index$All.Levins.std, HR.BD.Levins.index$All.Levins.std))

focal.compare.Levins.df$LS <- factor(focal.compare.Levins.df$LS,
                                levels = c("Konye - wet", "Konye - dry",
                                           "Ayos - wet", "Ayos - dry",
                                           "Bokito - wet", "Bokito - dry"))

y_expression <- expression("Levin's index" ~ "" ~(B[A]))

RA.compare.Levins.boxplot <- ggplot(focal.compare.Levins.df, aes(x = LS, y = Levins.index, fill = LS)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = y_expression) + ylim(0,0.6)

RA.compare.Levins.boxplot + annotate("text",
                                     x = 1:length(table(focal.compare.Levins.df$LS)),
                                     y = rep(0.57, times = 6),
                                     label = table(focal.compare.Levins.df$LS),
                                     col = "red",
                                     vjust = - 1)

focal.compare.Levins.df <- data.frame(LS = c(rep("HR - Konye - wet", times = length(HR.KW.Levins.index$All.Levins.std)),
                                              rep("HR - Konye - dry", times = length(HR.KD.Levins.index$All.Levins.std)),
                                              rep("HR - Ayos - wet", times = length(HR.AW.Levins.index$All.Levins.std)),
                                              rep("HR - Ayos - dry", times = length(HR.AD.Levins.index$All.Levins.std)),
                                              rep("HR - Bokito - wet", times = length(HR.BW.Levins.index$All.Levins.std)),
                                              rep("HR - Bokito - dry", times = length(HR.BD.Levins.index$All.Levins.std)),
                                              rep("RA - Konye - wet", times = length(RA.KW.Levins.index$All.Levins.std)),
                                              rep("RA - Konye - dry", times = length(RA.KD.Levins.index$All.Levins.std)),
                                              rep("RA - Ayos - wet", times = length(RA.AW.Levins.index$All.Levins.std)),
                                              rep("RA - Ayos - dry", times = length(RA.AD.Levins.index$All.Levins.std))),
                                       Counts = c(HR.KW.Levins.index$All.Levins.std, HR.KD.Levins.index$All.Levins.std, 
                                                  HR.AW.Levins.index$All.Levins.std, HR.AD.Levins.index$All.Levins.std, 
                                                  HR.BW.Levins.index$All.Levins.std, HR.BD.Levins.index$All.Levins.std,
                                                  RA.KW.Levins.index$All.Levins.std, RA.KD.Levins.index$All.Levins.std, 
                                                  RA.AW.Levins.index$All.Levins.std, RA.AD.Levins.index$All.Levins.std))

shapiro.test(focal.compare.Levins.df$Counts)
compare.Levins.kruskal <- kruskal.test(data = focal.compare.Levins.df, x = focal.compare.Levins.df$Counts,
                                        g = focal.compare.Levins.df$LS, formula = Counts ~ LS)

library(FSA)
dunnTest(Counts ~ LS, data = focal.compare.Levins.df, method = "bonferroni")


#####

# Pianka's niche overlap index
#####
library(EcoSimR)

pianka.null.function <- function(HR.df, RA.df) {
  null.model <- c(1:10000)
  
  for (n in null.model) {
    ra3.model <- ra3(speciesData = matrix(c(HR.df$wPOO, RA.df$wPOO), nrow = 2))
    
    null.model[n] <- pianka(m = ra3.model) 
  }
  return(null.model)
}

HR.KW.wPOO.pianka <- pianka(m = matrix(c(HR.KW.focal.df$wPOO, RA.KW.focal.df$wPOO), nrow = 2))
HR.KW.wPOO.null <- pianka.null.function(HR.KW.focal.df, RA.KW.focal.df)
wilcox.test(HR.KW.wPOO.null, mu = HR.KW.wPOO.pianka)

HR.KW.RRA.pianka <- pianka(m = matrix(c(HR.KW.focal.df$RRA, RA.KW.focal.df$RRA), nrow = 2))
HR.KW.RRA.null <- pianka.null.function(HR.KW.focal.df, RA.KW.focal.df)
wilcox.test(HR.KW.RRA.null, mu = HR.KW.RRA.pianka)

HR.KD.wPOO.pianka <- pianka(m = matrix(c(HR.KD.focal.df$wPOO, RA.KD.focal.df$wPOO), nrow = 2))
HR.KD.wPOO.null <- pianka.null.function(HR.KD.focal.df, RA.KD.focal.df)
wilcox.test(HR.KD.wPOO.null, mu = HR.KD.wPOO.pianka)

HR.KD.RRA.pianka <- pianka(m = matrix(c(HR.KD.focal.df$RRA, RA.KD.focal.df$RRA), nrow = 2))
HR.KD.RRA.null <- pianka.null.function(HR.KD.focal.df, RA.KD.focal.df)
wilcox.test(HR.KD.RRA.null, mu = HR.KD.RRA.pianka)

HR.AW.wPOO.pianka <- pianka(m = matrix(c(HR.AW.focal.df$wPOO, RA.AW.focal.df$wPOO), nrow = 2))
HR.AW.wPOO.null <- pianka.null.function(HR.AW.focal.df, RA.AW.focal.df)
wilcox.test(HR.AW.wPOO.null, mu = HR.AW.wPOO.pianka)

HR.AW.RRA.pianka <- pianka(m = matrix(c(HR.AW.focal.df$RRA, RA.AW.focal.df$RRA), nrow = 2))
HR.AW.RRA.null <- pianka.null.function(HR.AW.focal.df, RA.AW.focal.df)
wilcox.test(HR.AW.RRA.null, mu = HR.AW.RRA.pianka)

HR.AD.wPOO.pianka <- pianka(m = matrix(c(HR.AD.focal.df$wPOO, RA.AD.focal.df$wPOO), nrow = 2))
HR.AD.wPOO.null <- pianka.null.function(HR.AD.focal.df, RA.AD.focal.df)
wilcox.test(HR.AD.wPOO.null, mu = HR.AD.wPOO.pianka)

HR.AD.RRA.pianka <- pianka(m = matrix(c(HR.AD.focal.df$RRA, RA.AD.focal.df$RRA), nrow = 2))
HR.AD.RRA.null <- pianka.null.function(HR.AD.focal.df, RA.AD.focal.df)
wilcox.test(HR.AD.RRA.null, mu = HR.AD.RRA.pianka)

LocationSeasonList <- c("Konye - wet" , "Konye - dry" , "Ayos - wet" , "Ayos - dry")

pianka.data <- data.frame(LS = rep(LocationSeasonList, times = 2),
                          Metric = c(rep("wPOO", times = 4), rep("RRA", times = 4)),
                          Pianka = c(HR.KW.wPOO.pianka, HR.KD.wPOO.pianka,
                                     HR.AW.wPOO.pianka, HR.AD.wPOO.pianka,
                                     HR.KW.RRA.pianka, HR.KD.RRA.pianka,
                                     HR.AW.RRA.pianka, HR.AD.RRA.pianka)) 

pianka.data$LS <- factor(pianka.data$LS,
                                     levels = c("Konye - wet", "Konye - dry",
                                                "Ayos - wet", "Ayos - dry"))


pianka.plot <- ggplot(pianka.data, aes(x = LS, y = Pianka, fill = Metric)) + 
  geom_bar(stat="identity", position = "dodge") +
  theme_minimal() + scale_fill_manual(values = c("darkolivegreen3", "darkorange")) +
  theme(strip.background = element_blank(),
       strip.text = element_text(size = 12),
       axis.title.x = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.title.y = element_text(size = 12),
       axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5)) +
  labs(x = "Location - Season", y = "Pianka's index")

hist(pianka.data$Pianka)
shapiro.test(pianka.data$Pianka)
t.test(c(HR.KW.wPOO.pianka, HR.KD.wPOO.pianka,
                 HR.AW.wPOO.pianka, HR.AD.wPOO.pianka), 
            c(HR.KW.RRA.pianka, HR.KD.RRA.pianka,
                 HR.AW.RRA.pianka, HR.AD.RRA.pianka))

wilcox.test(c(HR.KW.wPOO.pianka, HR.KD.wPOO.pianka,
         HR.AW.wPOO.pianka, HR.AD.wPOO.pianka), 
       c(HR.KW.RRA.pianka, HR.KD.RRA.pianka,
         HR.AW.RRA.pianka, HR.AD.RRA.pianka))

wPOO.Pianka <- pianka.data[pianka.data$Metric=="RRA",]
mean(wPOO.Pianka$Pianka)

#####

# db-rda
# with presence/absence data
LS.new.zbj.focal <- new.zbj.focal
LS.new.zbj.focal <- cbind(LS.new.zbj.focal, paste(LS.new.zbj.focal$Location, "-", LS.new.zbj.focal$Season))
colnames(LS.new.zbj.focal) <- c(names(new.zbj.focal), "Location - Season")

capscale.function <- function(df) {
  
  focal.interact <- capscale(df[1:843] ~ df$Species, distance = "jaccard", add = TRUE)
  interact.output <- anova(focal.interact, step = 1000, perm.max = 200)
  
  vec = interact.output$SumOfSqs/sum(interact.output$SumOfSqs)*100
  table = interact.output
  table$SumOfSqs = vec
  table
  
  return(table)
}

KonyeWet.rda <- capscale.function(LS.new.zbj.focal[LS.new.zbj.focal$`Location - Season`=="Konye - wet",]) 
KonyeDry.rda <- capscale.function(LS.new.zbj.focal[LS.new.zbj.focal$`Location - Season`=="Konye - dry",]) 
AyosWet.rda <- capscale.function(LS.new.zbj.focal[LS.new.zbj.focal$`Location - Season`=="Ayos - wet",]) 
AyosDry.rda <- capscale.function(LS.new.zbj.focal[LS.new.zbj.focal$`Location - Season`=="Ayos - dry",]) 
BokitoWet.rda <- capscale.function(LS.new.zbj.focal[LS.new.zbj.focal$`Location - Season`=="Bokito - wet",]) 
BokitoDry.rda <- capscale.function(LS.new.zbj.focal[LS.new.zbj.focal$`Location - Season`=="Bokito - dry",]) 

# with RRA data
LS.new.zbj.rra <- new.zbj.rra
LS.new.zbj.rra <- cbind(LS.new.zbj.rra, paste(LS.new.zbj.rra$Location, "-", LS.new.zbj.rra$Season))
colnames(LS.new.zbj.rra) <- c(names(new.zbj.rra), "Location - Season")

capscale.function <- function(df) {
  
  focal.interact <- capscale(df[1:843] ~ df$Species, distance = "bray", add = TRUE)
  interact.output <- anova(focal.interact, step = 1000, perm.max = 200)
  
  vec = interact.output$SumOfSqs/sum(interact.output$SumOfSqs)*100
  table = interact.output
  table$SumOfSqs = vec
  table
  
  return(table)
}

KonyeWet.rda <- capscale.function(LS.new.zbj.rra[LS.new.zbj.rra$`Location - Season`=="Konye - wet",]) 
KonyeDry.rda <- capscale.function(LS.new.zbj.rra[LS.new.zbj.rra$`Location - Season`=="Konye - dry",]) 
AyosWet.rda <- capscale.function(LS.new.zbj.rra[LS.new.zbj.rra$`Location - Season`=="Ayos - wet",]) 
AyosDry.rda <- capscale.function(LS.new.zbj.rra[LS.new.zbj.rra$`Location - Season`=="Ayos - dry",]) 
BokitoWet.rda <- capscale.function(LS.new.zbj.rra[LS.new.zbj.rra$`Location - Season`=="Bokito - wet",]) 
BokitoDry.rda <- capscale.function(LS.new.zbj.rra[LS.new.zbj.rra$`Location - Season`=="Bokito - dry",]) 

# ordination plot for diet overlap between species
#####
library(vegan)                                                                              
library(dplyr)                                                                              
library(ggplot2)                                                                            
library(ggalt)                                                                              
library(ggforce)                                                                            
library(concaveman) 

# presence/absence NMDS
# between species
#####
new.zbj.focal <- new.zbj.focal[-c(which(new.zbj.focal$Location=="Bokito")),]
new.focal.info <- focal.info[-c(which(focal.info$Location=="Bokito")),]

focal.nmds <- metaMDS(new.zbj.focal[1:838], distance = "jaccard", k = 3, 
                      trymax=100, trace = FALSE)
plot(focal.nmds, main = paste("NMDS/Jaccard - Stress = ",                       
                              round(focal.nmds$stress, 3)))  
stressplot(focal.nmds)

focal.envfit <- envfit(focal.nmds, focal.info, permutations = 999, na.rm = TRUE)            
plot(focal.envfit)                                                                          

en_coord_cat <- as.data.frame(scores(focal.envfit, "factors")) 

data.scores <- as.data.frame(scores(focal.nmds))
data.scores <- cbind(data.scores, new.focal.info)
rownames(data.scores) <- c(1:152)

data.scores <- data.scores[-c(146, 149),]

focal.gg <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_mark_ellipse(expand=0, aes(fill=Species, colour = Species, size = Species), alpha = 0.1) +
  geom_point(data = data.scores, aes(shape = Species), size = 3, alpha = 0.5, show.legend = FALSE) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  scale_size_manual(values = c(0.2,0.2))

focal.gg
#####

# between each location-season
#####
new.data.scores <- as.data.frame(data.scores)
new.data.scores <- cbind(new.data.scores, paste(new.data.scores$Location, "-", new.data.scores$Season))
colnames(new.data.scores) <- c("NMDS1", "NMDS2", "NMDS3", "Species", "Location", "Season", "Location - Season")

new.data.scores$`Location - Season` <- factor(new.data.scores$`Location - Season`,
                                              levels = c("Konye - wet", "Konye - dry",
                                                         "Ayos - wet", "Ayos - dry"))

ordination.label <- c()
for (LS in LS.list) {
  NDS <- new.data.scores[new.data.scores$`Location - Season`==LS,]
  
  for (sp in unique(new.data.scores$Species)) {
    NDS.species <- NDS[NDS$Species==sp,]
    ordination.label <- append(ordination.label, length(NDS.species$Species))
  }
}

label.df <- data.frame(Species = rep(unique(new.data.scores$Species), times = 3),
                       Place = LS.list,
                       Count = ordination.label)

focal.gg2 <- ggplot(data = new.data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_mark_ellipse(expand=0, aes(fill=Species, colour = Species, size = Species), alpha = 0.1) +
  geom_point(data = new.data.scores, aes(shape = Species,), size = 3, alpha = 0.5, show.legend = FALSE) + 
  scale_size_manual(values = c(0.2,0.2)) + facet_wrap(vars(`Location - Season`), ncol = 2) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"),
        panel.spacing = unit(0.5, "cm", data = NULL),
        strip.background = element_blank()) 

focal.gg2
#####

# NMDS with wPOO summaries
#####
wPOO.matrix <- matrix(c(HR.KW.focal.df$wPOO, HR.KD.focal.df$wPOO,
                        HR.AW.focal.df$wPOO, HR.AD.focal.df$wPOO,
                        RA.KW.focal.df$wPOO, RA.KD.focal.df$wPOO,
                        RA.AW.focal.df$wPOO, RA.AD.focal.df$wPOO), nrow = 8)

wPOO.info <- data.frame(Species = c(rep("H. ruber", times = 4),
                                    rep("R. alcyone", times = 4)),
                        LS = rep(c("Konye - wet", "Konye - dry",
                               "Ayos - wet", "Ayos - dry"), times = 2))

wPOO.nmds <- metaMDS(wPOO.matrix, distance = "jaccard", k = 2, 
                      trymax=100, trace = FALSE)

plot(wPOO.nmds, main = paste("NMDS/Jaccard - Stress = ",                       
                              round(wPOO.nmds$stress, 3)))  
stressplot(wPOO.nmds)

wPOO.envfit <- envfit(wPOO.nmds, wPOO.info, permutations = 999, na.rm = TRUE)            
plot(wPOO.envfit)                                                                          

en_coord_cat <- as.data.frame(scores(focal.envfit, "factors"))

wPOO.data.scores <- as.data.frame(scores(wPOO.nmds))
wPOO.data.scores <- cbind(data.scores, wPOO.info) 

wPOO.nmds.graph <- ggplot(data = wPOO.data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_mark_ellipse(expand=0, aes(fill=Species, colour = Species), alpha = 0.1) +
  geom_point(data = wPOO.data.scores, aes(shape = Species), size = 3, alpha = 0.5, show.legend = FALSE) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) 

wPOO.nmds.graph
#####

# RRA NMDS
# between species
#####
new.zbj.rra <- new.zbj.rra[-c(which(new.zbj.rra$Location=="Bokito")),]

focal.nmds <- metaMDS(new.zbj.rra[1:838], distance = "jaccard", k = 3, 
                      trymax=100, trace = FALSE)
plot(focal.nmds, main = paste("NMDS/Jaccard - Stress = ",                       
                              round(focal.nmds$stress, 3)))  
stressplot(focal.nmds)

focal.envfit <- envfit(focal.nmds, new.focal.info, permutations = 999, na.rm = TRUE)            
plot(focal.envfit)                                                                          

en_coord_cat <- as.data.frame(scores(focal.envfit, "factors")) 

data.scores <- as.data.frame(scores(focal.nmds))
data.scores <- cbind(data.scores, new.focal.info)                                               
rownames(data.scores) <- c(1:152)

data.scores <- data.scores[-c(146, 149),]


focal.gg <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_mark_ellipse(expand=0, aes(fill=Species, colour = Species, size = Species), alpha = 0.1) +
  geom_point(data = data.scores, aes(shape = Species), size = 3, alpha = 0.5, show.legend = FALSE) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  scale_size_manual(values = c(0.2,0.2))

focal.gg
#####

# between each location-season
#####
new.data.scores <- as.data.frame(data.scores)
new.data.scores <- cbind(new.data.scores, paste(new.data.scores$Location, "-", new.data.scores$Season))
colnames(new.data.scores) <- c("NMDS1", "NMDS2", "NMDS3", "Species", "Location", "Season", "Location - Season")

new.data.scores$`Location - Season` <- factor(new.data.scores$`Location - Season`,
                                              levels = c("Konye - wet", "Konye - dry",
                                                         "Ayos - wet", "Ayos - dry",
                                                         "Bokito - wet", "Bokito - dry"))

ordination.label <- c()
for (LS in LS.list) {
  NDS <- new.data.scores[new.data.scores$`Location - Season`==LS,]
  
  for (sp in unique(new.data.scores$Species)) {
    NDS.species <- NDS[NDS$Species==sp,]
    ordination.label <- append(ordination.label, length(NDS.species$Species))
  }
}

label.df <- data.frame(Species = rep(unique(new.data.scores$Species), times = 3),
                       Place = LS.list,
                       Count = ordination.label)

focal.gg2 <- ggplot(data = new.data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_mark_ellipse(expand=0, aes(fill=Species, colour = Species, size = Species), alpha = 0.1) +
  geom_point(data = new.data.scores, aes(shape = Species,), size = 3, alpha = 0.5, show.legend = FALSE) + 
  scale_size_manual(values = c(0.2,0.2)) + facet_wrap(vars(`Location - Season`), ncol = 2) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"),
        panel.spacing = unit(0.5, "cm", data = NULL),
        strip.background = element_blank()) 

focal.gg2
#####

# NMDS with RRA summaries
#####
RRA.matrix <- matrix(c(HR.KW.focal.df$RRA, HR.KD.focal.df$RRA,
                        HR.AW.focal.df$RRA, HR.AD.focal.df$RRA,
                        RA.KW.focal.df$RRA, RA.KD.focal.df$RRA,
                        RA.AW.focal.df$RRA, RA.AD.focal.df$RRA), nrow = 8)

RRA.info <- data.frame(Species = c(rep("H. ruber", times = 4),
                                    rep("R. alcyone", times = 4)),
                        LS = rep(c("Konye - wet", "Konye - dry",
                                   "Ayos - wet", "Ayos - dry"), times = 2))

RRA.nmds <- metaMDS(RRA.matrix, distance = "jaccard", k = 2, 
                     trymax=100, trace = FALSE)

plot(RRA.nmds, main = paste("NMDS/Jaccard - Stress = ",                       
                             round(RRA.nmds$stress, 3)))  
stressplot(RRA.nmds)

RRA.data.scores <- as.data.frame(scores(RRA.nmds))
RRA.data.scores <- cbind(data.scores, RRA.info) 

RRA.nmds.graph <- ggplot(data = RRA.data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_mark_ellipse(expand=0, aes(fill=Species, colour = Species), alpha = 0.1) +
  geom_point(data = RRA.data.scores, aes(shape = Species), size = 3, alpha = 0.5, show.legend = FALSE) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) 

RRA.nmds.graph
#####

# Occurrence of pest sequences
pest.data <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\Pest sequences.xlsx")
pest.data <- as.data.frame(pest.data)
pest.data <- pest.data[,-c(56,185,272,274:275, 277:291)]

rownames(pest.data) <- c(pest.data$OTU)

pest.data <- cbind(pest.data, paste(pest.data$Pest, "-", pest.data$sacc))

#pest.OTUs <- pest.data[c(6,13,26,27,57,62,73,80:82,91,97:102,
#                         111:114,121,129:137,142,149,170:183,
#                         197:202,207:209,15,74,75,88,107,122,
#                         150,203,210,211), "OTU"]

#pest.nbrs <- c(which(pest.data$OTU %in% pest.OTUs))

HR.pest.df <- pest.data[,c(c(HR.nbrs), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]

RA.pest.df <- pest.data[,c(c(RA.nbrs), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]

HR.KW.pest.df <- HR.pest.df[,c(as.character(c(HR.KW.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
HR.KW.pest.df <- cbind(HR.KW.pest.df, data.frame(HR.KW.focal.df[c(which(HR.KW.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
HR.KW.pest.df.cut <- filter.rows(HR.KW.pest.df[,1:13], 1, z)
HR.KW.pest.df <- HR.KW.pest.df[-c(HR.KW.pest.df.cut),]

HR.KD.pest.df <- HR.pest.df[,c(as.character(c(HR.KD.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
HR.KD.pest.df <- cbind(HR.KD.pest.df, data.frame(HR.KD.focal.df[c(which(HR.KD.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
HR.KD.pest.df.cut <- filter.rows(HR.KD.pest.df[,1:24], 1, z)
HR.KD.pest.df <- HR.KD.pest.df[-c(HR.KD.pest.df.cut),]

HR.AW.pest.df <- HR.pest.df[,c(as.character(c(HR.AW.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
HR.AW.pest.df <- cbind(HR.AW.pest.df, data.frame(HR.AW.focal.df[c(which(HR.AW.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
HR.AW.pest.df.cut <- filter.rows(HR.AW.pest.df[,1:10], 1, z)
HR.AW.pest.df <- HR.AW.pest.df[-c(HR.AW.pest.df.cut),]

HR.AD.pest.df <- HR.pest.df[,c(as.character(c(HR.AD.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
HR.AD.pest.df <- cbind(HR.AD.pest.df, data.frame(HR.AD.focal.df[c(which(HR.AD.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
HR.AD.pest.df.cut <- filter.rows(HR.AD.pest.df[,1:25], 1, z)
HR.AD.pest.df <- HR.AD.pest.df[-c(HR.AD.pest.df.cut),]

HR.BW.pest.df <- HR.pest.df[,c(as.character(c(HR.BW.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
HR.BW.pest.df <- cbind(HR.BW.pest.df, data.frame(HR.BW.focal.df[c(which(HR.BW.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
HR.BW.pest.df.cut <- filter.rows(HR.BW.pest.df[,1:4], 1, z)
HR.BW.pest.df <- HR.BW.pest.df[-c(HR.BW.pest.df.cut),]

HR.BD.pest.df <- HR.pest.df[,c(as.character(c(HR.BD.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
HR.BD.pest.df <- cbind(HR.BD.pest.df, data.frame(HR.BD.focal.df[c(which(HR.BD.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
HR.BD.pest.df.cut <- filter.rows(HR.BD.pest.df[,1:3], 1, z)
HR.BD.pest.df <- HR.BD.pest.df[-c(HR.BD.pest.df.cut),]

RA.KW.pest.df <- RA.pest.df[,c(as.character(c(RA.KW.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
RA.KW.pest.df <- cbind(RA.KW.pest.df, data.frame(RA.KW.focal.df[c(which(RA.KW.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
RA.KW.pest.df.cut <- filter.rows(RA.KW.pest.df[,1:29], 1, z)
RA.KW.pest.df <- RA.KW.pest.df[-c(RA.KW.pest.df.cut),]

RA.KD.pest.df <- RA.pest.df[,c(as.character(c(RA.KD.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
RA.KD.pest.df <- cbind(RA.KD.pest.df, data.frame(RA.KD.focal.df[c(which(RA.KD.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
RA.KD.pest.df.cut <- filter.rows(RA.KD.pest.df[,1:18], 1, z)
RA.KD.pest.df <- RA.KD.pest.df[-c(RA.KD.pest.df.cut),]

RA.AW.pest.df <- RA.pest.df[,c(as.character(c(RA.AW.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
RA.AW.pest.df <- cbind(RA.AW.pest.df, data.frame(RA.AW.focal.df[c(which(RA.AW.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
RA.AW.pest.df.cut <- filter.rows(RA.AW.pest.df[,1:20], 1, z)
RA.AW.pest.df <- RA.AW.pest.df[-c(RA.AW.pest.df.cut),]

RA.AD.pest.df <- RA.pest.df[,c(as.character(c(RA.AD.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
RA.AD.pest.df <- cbind(RA.AD.pest.df, data.frame(RA.AD.focal.df[c(which(RA.AD.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
RA.AD.pest.df.cut <- filter.rows(RA.AD.pest.df[,1:13], 1, z)
RA.AD.pest.df <- RA.AD.pest.df[-c(RA.AD.pest.df.cut),]

RA.BW.pest.df <- RA.pest.df[,c(as.character(c(RA.BW.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
RA.BW.pest.df <- cbind(RA.BW.pest.df, data.frame(RA.BW.focal.df[c(which(RA.BW.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
RA.BW.pest.df.cut <- which(RA.BW.pest.df[,1]==0)
RA.BW.pest.df <- RA.BW.pest.df[-c(RA.BW.pest.df.cut),]

RA.BD.pest.df <- RA.pest.df[,c(as.character(c(RA.BD.nbrs)), "paste(pest.data$Pest, \"-\", pest.data$sacc)", "order", "family", "genus", "species")]
RA.BD.pest.df <- cbind(RA.BD.pest.df, data.frame(RA.BD.focal.df[c(which(RA.BD.focal.df$OTU %in% pest.data$OTU)),
                                                                c("wPOO", "RRA")]))
RA.BD.pest.df.cut <- which(RA.BD.pest.df[,1]==0)
RA.BD.pest.df <- RA.BD.pest.df[-c(RA.BD.pest.df.cut),]

saccList <- unique(pest.data$`paste(pest.data$Pest, \"-\", pest.data$sacc)`)
# sacc summary function
#####
saccSummary <- function(df) {
  wPOOsaccMean <- c()
  RRAsaccMean <- c()
  
  for (sacc in saccList) {
    sacc.df <- df[df$`paste(pest.data$Pest, \"-\", pest.data$sacc)`==sacc,]
    
    ifelse(nrow(sacc.df) > 0,
           wPOOsaccMean <- append(wPOOsaccMean, mean(sacc.df$wPOO)),
           wPOOsaccMean <- append(wPOOsaccMean, 0))
    
    ifelse(nrow(sacc.df) > 0,
           RRAsaccMean <- append(RRAsaccMean, mean(sacc.df$RRA)),
           RRAsaccMean <- append(RRAsaccMean, 0))
  }
  return(list(wPOOsaccMean, RRAsaccMean))
}
#####

HR.KW.sacc <- saccSummary(HR.KW.pest.df)
HR.KD.sacc <- saccSummary(HR.KD.pest.df)
HR.AW.sacc <- saccSummary(HR.AW.pest.df)
HR.AD.sacc <- saccSummary(HR.AD.pest.df)
HR.BW.sacc <- saccSummary(HR.BW.pest.df)
HR.BD.sacc <- saccSummary(HR.BD.pest.df)

RA.KW.sacc <- saccSummary(RA.KW.pest.df)
RA.KD.sacc <- saccSummary(RA.KD.pest.df)
RA.AW.sacc <- saccSummary(RA.AW.pest.df)
RA.AD.sacc <- saccSummary(RA.AD.pest.df)
RA.BW.sacc <- list(rep(0, times = length(saccList)),
                   rep(0, times = length(saccList)))
RA.BD.sacc <- list(rep(0, times = length(saccList)),
                   rep(0, times = length(saccList)))


saccSummary.df <- data.frame(Species = rep(c("H.ruber", "R.alcyone"), each = 6*length(saccList), times = 2),
                             LS = rep(c("Konye - wet", "Konye - dry", "Ayos - wet",
                                        "Ayos - dry", "Bokito - wet", "Bokito - dry"),
                                      each = length(saccList), times = 4),
                             Sacc = rep(saccList, times = 24),
                             Summary = rep(c("RRA", "wPOO"), each = 12*length(saccList)),
                             Value = c(HR.KW.sacc[[2]], HR.KD.sacc[[2]], HR.AW.sacc[[2]],
                                       HR.AD.sacc[[2]], HR.BW.sacc[[2]], HR.BD.sacc[[2]],
                                       RA.KW.sacc[[2]], RA.KD.sacc[[2]], RA.AW.sacc[[2]],
                                       RA.AD.sacc[[2]], RA.BW.sacc[[2]], RA.BD.sacc[[2]],
                                       HR.KW.sacc[[1]], HR.KD.sacc[[1]], HR.AW.sacc[[1]],
                                      HR.AD.sacc[[1]], HR.BW.sacc[[1]], HR.BD.sacc[[1]],
                                      RA.KW.sacc[[1]], RA.KD.sacc[[1]], RA.AW.sacc[[1]],
                                      RA.AD.sacc[[1]], RA.BW.sacc[[1]], RA.BD.sacc[[1]]))

saccSummary.cut <- which(saccSummary.df$Value < 0.00001)
saccSummary.df <- saccSummary.df[-c(saccSummary.cut),]

pestPlot <- ggplot(saccSummary.df[saccSummary.df$LS=="Bokito - dry",], aes(x = Sacc, y = Value, 
                                                                             fill = factor(Summary, levels = c("wPOO", "RRA")))) +
  geom_bar(stat = "identity", position = "stack", ) + theme_minimal() +
  facet_wrap(vars(Species), ncol = 1) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  labs(x = "Pest", y = "Average abundance", fill = "Metric") +
  scale_fill_manual(values = c("darkorange", "darkolivegreen3"))

RApestPlot <- ggplot(saccSummary.df[saccSummary.df$Species=="R.alcyone",], aes(x = Sacc, y = Value, fill = Summary)) +
  geom_bar(stat = "identity", position = "stack", ) + theme_minimal() +
  facet_wrap(vars(factor(LS, levels = c("Konye - wet", "Konye - dry", "Ayos - wet",
                                        "Ayos - dry", "Bokito - wet", "Bokito - dry")))) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 1),
        axis.text.y = element_text(size = 8)) +
  labs(x = "Pest", y = "Average abundance", fill = "Metric")


# Venn diagram comparing pest consumption between species
#####
library(RAM)
compare.pest.OTUs <- list(H.ruber = HR.pest.df$OTU,
                          R.alcyone = RA.pest.df$OTU)

group.venn(compare.pest.OTUs, label=FALSE,  lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(330, 30), cat.dist = 0.05)

#####

# Compare average number of pest sequences between species
#####
HR.pest.OTUs <- count.OTUs(HR.pest.df[,1:79])
RA.pest.OTUs <- count.OTUs(RA.pest.df[,1:82])

HR.pest.df <- HR.pest.df[,-c(which(HR.pest.OTUs==0))]
RA.pest.df <- RA.pest.df[,-c(which(RA.pest.OTUs==0))]

HR.pest.OTUs <- count.OTUs(HR.pest.df[,1:66])
RA.pest.OTUs <- count.OTUs(RA.pest.df[,1:65])

compare.pest.OTUs.df <- data.frame(Species = c(rep("H.ruber", times = length(HR.pest.OTUs)),
                                               rep("R.alcyone", times = length(RA.pest.OTUs))),
                                   Counts = c(HR.pest.OTUs, RA.pest.OTUs))

focal.compare.OTUs.boxplot <- ggplot(compare.pest.OTUs.df, aes(x = Species, y = Counts, fill = Species)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.text.y = element_text(size = 10)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Species", y = "Number of OTUs")

shapiro.test(compare.pest.OTUs.df$Counts)
wilcox.test(HR.pest.OTUs, RA.pest.OTUs)
#####

HR.KW.pest.OTUs <- count.OTUs(HR.KW.pest.df[,1:13])
HR.KD.pest.OTUs <- count.OTUs(HR.KD.pest.df[,1:24])
HR.AW.pest.OTUs <- count.OTUs(HR.AW.pest.df[,1:10])
HR.AD.pest.OTUs <- count.OTUs(HR.AD.pest.df[,1:25])
HR.BW.pest.OTUs <- count.OTUs(HR.BW.pest.df[,1:4])
HR.BD.pest.OTUs <- count.OTUs(HR.BD.pest.df[,1:3])

RA.KW.pest.OTUs <- count.OTUs(RA.KW.pest.df[,1:29])
RA.KD.pest.OTUs <- count.OTUs(RA.KD.pest.df[,1:18])
RA.AW.pest.OTUs <- count.OTUs(RA.AW.pest.df[,1:20])
RA.AD.pest.OTUs <- count.OTUs(RA.AD.pest.df[,1:13])
RA.BW.pest.OTUs <- length(which(RA.BW.pest.df[,1]==1))
RA.BD.pest.OTUs <- length(which(RA.BD.pest.df[,1]==1))

# How does number of pest average number of pest sequences consumed vary between locations and seasons
#####
HR.compare.pest.OTUs.df <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.KW.pest.OTUs, HR.KD.pest.OTUs,
                                                                                  HR.AW.pest.OTUs, HR.AD.pest.OTUs,
                                                                                  HR.BW.pest.OTUs, HR.BD.pest.OTUs)))),
                                      LS = c(rep("Konye - wet", times = length(HR.KW.pest.OTUs)), 
                                             rep("Konye - dry", times = length(HR.KD.pest.OTUs)),
                                             rep("Ayos - wet", times = length(HR.AW.pest.OTUs)), 
                                             rep("Ayos - dry", times = length(HR.AD.pest.OTUs)),
                                             rep("Bokito - wet", times = length(HR.BW.pest.OTUs)), 
                                             rep("Bokito - dry", times = length(HR.BD.pest.OTUs))),
                                      Counts = c(HR.KW.pest.OTUs, HR.KD.pest.OTUs,
                                                 HR.AW.pest.OTUs, HR.AD.pest.OTUs,
                                                 HR.BW.pest.OTUs, HR.BD.pest.OTUs))

HR.compare.pest.OTUs.df.boxplot <- ggplot(HR.compare.pest.OTUs.df, aes(x = LS, y = Counts, fill = LS)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 90)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = "Number of pest OTUs") 

HR.compare.pest.OTUs.df.boxplot + annotate("text",
                                           x = 1:length(table(HR.compare.pest.OTUs.df$LS)),
                                           y = rep(4, times = 6),
                                           label = table(HR.compare.pest.OTUs.df$LS),
                                           col = "red",
                                           vjust = - 1)

library(FSA)
HR.compare.pest.OTUs.kruskal <- kruskal.test(data = HR.compare.pest.OTUs.df, x = HR.compare.pest.OTUs.df$Counts,
                                             g = HR.compare.pest.OTUs.df$LS, formula = Counts ~ LS)

dunnTest(Counts ~ LS, data = HR.compare.pest.OTUs.df, method = "bonferroni")

RA.compare.pest.OTUs.df <- data.frame(Species = c(rep("R.alcyone", times = length(c(RA.KW.pest.OTUs, RA.KD.pest.OTUs,
                                                                                    RA.AW.pest.OTUs, RA.AD.pest.OTUs,
                                                                                    RA.BW.pest.OTUs, RA.BD.pest.OTUs)))),
                                      LS = c(rep("Konye - wet", times = length(RA.KW.pest.OTUs)), 
                                             rep("Konye - dry", times = length(RA.KD.pest.OTUs)),
                                             rep("Ayos - wet", times = length(RA.AW.pest.OTUs)), 
                                             rep("Ayos - dry", times = length(RA.AD.pest.OTUs)),
                                             rep("Bokito - wet", times = length(RA.BW.pest.OTUs)), 
                                             rep("Bokito - dry", times = length(RA.BD.pest.OTUs))),
                                      Counts = c(RA.KW.pest.OTUs, RA.KD.pest.OTUs,
                                                 RA.AW.pest.OTUs, RA.AD.pest.OTUs,
                                                 RA.BW.pest.OTUs, RA.BD.pest.OTUs))

RA.compare.pest.OTUs.df.boxplot <- ggplot(RA.compare.pest.OTUs.df, aes(x = LS, y = Counts, fill = LS)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10, angle = 90)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location - Season", y = "Number of pest OTUs") 

RA.compare.pest.OTUs.df.boxplot + annotate("text",
                                           x = 1:length(table(RA.compare.pest.OTUs.df$LS)),
                                           y = rep(5, times = 6),
                                           label = table(RA.compare.pest.OTUs.df$LS),
                                           col = "red",
                                           vjust = - 1)

RA.compare.pest.OTUs.kruskal <- kruskal.test(data = RA.compare.pest.OTUs.df, x = RA.compare.pest.OTUs.df$Counts,
                                             g = RA.compare.pest.OTUs.df$LS, formula = Counts ~ LS)

dunnTest(Counts ~ LS, data = RA.compare.pest.OTUs.df, method = "bonferroni")
#####

# Occurrence of cocoa pest sequences
#####
find.pests <- function(df) {
  pest.rows <- c()
  for (row in 1:nrow(df)) {
    if (df[row, "OTU"] %in% pest.OTUs) {
      pest.rows <- append(pest.rows, row)
    }
  } 
  return(pest.rows)
}

HR.KW.cocoa.pests <- find.pests(HR.KW.pest.df)
HR.KW.cocoa.pests <- HR.KW.pest.df[c(HR.KW.cocoa.pests),]

HR.KD.cocoa.pests <- find.pests(HR.KD.pest.df)
HR.KD.cocoa.pests <- HR.KD.pest.df[c(HR.KD.cocoa.pests),]

HR.AW.cocoa.pests <- find.pests(HR.AW.pest.df)
HR.AW.cocoa.pests <- HR.AW.pest.df[c(HR.AW.cocoa.pests),]

HR.AD.cocoa.pests <- find.pests(HR.AD.pest.df)
HR.AD.cocoa.pests <- HR.AD.pest.df[c(HR.AD.cocoa.pests),]

##

HR.KW.cut <- filter.rows(HR.KW.pest.df[,1:13], 1, z)
HR.KW.pest.df <- HR.KW.pest.df[-c(HR.KW.cut),]

HR.KD.cut <- filter.rows(HR.KD.pest.df[,1:24], 1, z)
HR.KD.pest.df <- HR.KD.pest.df[-c(HR.KD.cut),]

HR.AW.cut <- filter.rows(HR.AW.pest.df[,1:10], 1, z)
HR.AW.pest.df <- HR.AW.pest.df[-c(HR.AW.cut),]

HR.AD.cut <- filter.rows(HR.AD.pest.df[,1:25], 1, z)
HR.AD.pest.df <- HR.AD.pest.df[-c(HR.AD.cut),]

HR.BW.cut <- filter.rows(HR.BW.pest.df[,1:4], 1, z)
HR.BW.pest.df <- HR.BW.pest.df[-c(HR.BW.cut),]

HR.BD.cut <- filter.rows(HR.BD.pest.df[,1:3], 1, z)
HR.BD.pest.df <- HR.BD.pest.df[-c(HR.BD.cut),]


RA.KW.cut <- filter.rows(RA.KW.pest.df[,1:29], 1, z)
RA.KW.pest.df <- RA.KW.pest.df[-c(RA.KW.cut),]

RA.KD.cut <- filter.rows(RA.KD.pest.df[,1:18], 1, z)
RA.KD.pest.df <- RA.KD.pest.df[-c(RA.KD.cut),]

RA.AW.cut <- filter.rows(RA.AW.pest.df[,1:20], 1, z)
RA.AW.pest.df <- RA.AW.pest.df[-c(RA.AW.cut),]

RA.AD.cut <- filter.rows(RA.AD.pest.df[,1:13], 1, z)
RA.AD.pest.df <- RA.AD.pest.df[-c(RA.AD.cut),]

RA.BD.pest.df <- RA.BD.pest.df[c(which(RA.BD.pest.df[,1]==1)),]
#####

# Occurrence of cocoa pest sequences
find.pests <- function(df) {
  pest.rows <- c()
  for (row in 1:nrow(df)) {
    if (df[row, "OTU"] %in% pest.OTUs) {
      pest.rows <- append(pest.rows, row)
    }
  } 
  return(pest.rows)
}

HR.KW.cocoa.pests <- find.pests(HR.KW.pest.df)
HR.KW.cocoa.pests <- HR.KW.pest.df[c(HR.KW.cocoa.pests),]

HR.KD.cocoa.pests <- find.pests(HR.KD.pest.df)
HR.KD.cocoa.pests <- HR.KD.pest.df[c(HR.KD.cocoa.pests),]

HR.AW.cocoa.pests <- find.pests(HR.AW.pest.df)
HR.AW.cocoa.pests <- HR.AW.pest.df[c(HR.AW.cocoa.pests),]

HR.AD.cocoa.pests <- find.pests(HR.AD.pest.df)
HR.AD.cocoa.pests <- HR.AD.pest.df[c(HR.AD.cocoa.pests),]

######

# is there a difference in the time the bats were caught
light.data <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\Copy of Light conditions.xlsx")
light.data <- as.data.frame(light.data)

HR.light <- light.data[c(which(light.data$Species=="Hrub")),]
HR.Konye.light <- HR.light[HR.light$Location=="Konye",]
HR.Ayos.light <- HR.light[HR.light$Location=="Ayos",]
HR.Bokito.light <- HR.light[HR.light$Location=="Bokito",]

RA.light <- light.data[c(which(light.data$Species=="Ralc")),]
RA.Konye.light <- RA.light[RA.light$Location=="Konye",]
RA.Ayos.light <- RA.light[RA.light$Location=="Ayos",]
RA.Bokito.light <- RA.light[RA.light$Location=="Bokito",]

nightTime <- c("18hr", "19hr", "20hr", "21hr", "22hr", "23hr", "00hr", "1hr")

lightFreq.function <- function(df) {
  timeFreq <- c()
  for (time in nightTime) {
    timeFreq <- append(timeFreq, length(which(df$Time==time)))
  }
  return(timeFreq)
}

allTimeFreq <- data.frame(Species = rep(c("H. ruber", "R. alcyone"), each = length(nightTime)*3),
                          Location = rep(c("Konye", "Ayos", "Bokito"), each = length(nightTime), times = 2),
                          Time = rep(nightTime, times = 6),
                          Freq. = c(lightFreq.function(HR.Konye.light),
                                    lightFreq.function(HR.Ayos.light),
                                    lightFreq.function(HR.Bokito.light),
                                    lightFreq.function(RA.Konye.light),
                                    lightFreq.function(RA.Ayos.light),
                                    lightFreq.function(RA.Bokito.light)))

allTimeFreq$Time <- factor(allTimeFreq$Time, levels = nightTime)

light.plot <- ggplot(allTimeFreq, aes(x = Time, y = Freq., fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") + theme_minimal()  +
  facet_wrap(vars(Location))



focal.light.data <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\Light conditions.xlsx")
focal.light.data <- as.data.frame(focal.light.data)


# diet composition 
# function to find most common Lepidopteran families in each location and season
findLepidoptera <- function(df) {
  lepidoptera_df <- df[df$order=="Lepidoptera",]
  lepidoptera_fam <- unique(lepidoptera_df$family)
  fam_na <- which(lepidoptera_fam=="NA")
  
  return(lepidoptera_fam[-fam_na])
}
#####

RA.KW.focal.cut <- cbind(RA.KW.focal.cut, insect.data.order[-c(RA.KW.focal.trim),])
colnames(RA.KW.focal.cut) <- c(RA.KW.nbrs, "OTU", "order", "family", "genus", "species")

RA.KD.focal.cut <- cbind(RA.KD.focal.cut, insect.data.order[-c(RA.KD.focal.trim),])
colnames(RA.KD.focal.cut) <- c(RA.KD.nbrs, "OTU", "order", "family", "genus", "species")

RA.AW.focal.cut <- cbind(RA.AW.focal.cut, insect.data.order[-c(RA.AW.focal.trim),])
colnames(RA.AW.focal.cut) <- c(RA.AW.nbrs, "OTU", "order", "family", "genus", "species")

RA.AD.focal.cut <- cbind(RA.AD.focal.cut, insect.data.order[-c(RA.AD.focal.trim),])
colnames(RA.AD.focal.cut) <- c(RA.AD.nbrs, "OTU", "order", "family", "genus", "species")

RA.BW.focal.cut <- cbind(RA.BW.focal.cut, insect.data.order[-c(RA.BW.focal.trim),])
colnames(RA.BW.focal.cut) <- c(RA.BW.nbrs, "OTU", "order", "family", "genus", "species")

RA.BD.focal.cut <- cbind(RA.BD.focal.cut, insect.data.order[-c(RA.BD.focal.trim),])
colnames(RA.BD.focal.cut) <- c(RA.BD.nbrs, "OTU", "order", "family", "genus", "species")

RA.KW.fam <- findLepidoptera(RA.KW.focal.cut)
RA.KD.fam <- findLepidoptera(RA.KD.focal.cut)
RA.AW.fam <- findLepidoptera(RA.AW.focal.cut)
RA.AD.fam <- findLepidoptera(RA.AD.focal.cut)
RA.BW.fam <- findLepidoptera(RA.BW.focal.cut)
RA.BD.fam <- findLepidoptera(RA.BD.focal.cut)

RA_Lepidoptera <- list("Konye-wet" = RA.KW.fam,
                       "Konye-dry" = RA.KD.fam,
                       "Ayos-wet" = RA.AW.fam,
                       "Ayos-dry" = RA.AD.fam,
                       "Bokito-wet" = RA.BW.fam,
                       "Bokito-dry" = RA.BD.fam)

allLepidoptera <- insect.data.order[insect.data.order$order=="Lepidoptera",]
allLepidoptera <- unique(allLepidoptera$family)
allLepidoptera <- allLepidoptera[-c(which(allLepidoptera=="NA"))]

lepidopteraYesNo <- matrix(ncol = 0, nrow = 35)

for (l in RA_Lepidoptera) {
  YesNo <- c()
  for (all in allLepidoptera) {
    ifelse(all %in% l, YesNo <- append(YesNo, 1), YesNo <- append(YesNo, 0))
  }
  lepidopteraYesNo <- cbind(lepidopteraYesNo, YesNo)
}

colnames(lepidopteraYesNo) <- c("Konye-wet", "Konye-dry", "Ayos-wet",
                                "Ayos-dry", "Bokito-wet", "Bokito-dry")
rownames(lepidopteraYesNo) <- allLepidoptera

lepidopteraCut <- filter.rows(lepidopteraYesNo, 1, z)
lepidopteraYesNo <- lepidopteraYesNo[-c(lepidopteraCut),]

rowSum <- c()
for (row in 1:nrow(lepidopteraYesNo)) {
  rowSum <- append(rowSum, sum(lepidopteraYesNo[row,]))
}

lepidopteraYesNo <- cbind(lepidopteraYesNo, rowSum)
#lepidopteraYesNo <- lepidopteraYesNo[order(lepidopteraYesNo[,7]),]
#lepidopteraYesNo <- lepidopteraYesNo[,-7]

allLepidopteraCut <- allLepidoptera[-c(lepidopteraCut)]

lepidopteraYesNo <- cbind(lepidopteraYesNo, allLepidopteraCut)

lepidopteraYesNo <- lepidopteraYesNo[order(factor(lepidopteraYesNo[,8], 
                                                   levels = c("Noctuidae","Pyralidae", "Geometridae", "Tortricidae",
                                                              "Tineidae", "Crambidae", "Lecithoceridae",
                                                              "Oecophoridae", "Cosmopterigidae", "Erebidae", "Gelechiidae",
                                                              "Limacodidae", "Nolidae", "Autostichidae",
                                                              "Anthelidae", "Carposinidae", "Blastobasidae",
                                                              "Yponomeutidae", "Cossidae",
                                                              "Saturniidae", "Choreutidae",
                                                              "Scythrididae", "Lymantriidae",
                                                              "Coleophoridae", "Notodontidae", "Lasiocampidae",
                                                              "Eupterotidae","Elachistidae", "Euteliidae"))),]

lepidopteraYesNo <- lepidopteraYesNo[,-c(7,8)]
new.lepidopteraYesNo <- matrix(ncol = 0, nrow = 29)

for (col in 1:ncol(lepidopteraYesNo)) {
  new.lepidopteraYesNo <- cbind(new.lepidopteraYesNo, as.numeric(lepidopteraYesNo[,col]))
}

rownames(new.lepidopteraYesNo) <- c("Noctuidae","Pyralidae", "Geometridae", "Tortricidae",
                                    "Tineidae", "Crambidae", "Lecithoceridae",
                                    "Oecophoridae", "Cosmopterigidae", "Erebidae", "Gelechiidae",
                                    "Limacodidae", "Nolidae", "Autostichidae",
                                    "Anthelidae", "Carposinidae", "Blastobasidae",
                                    "Yponomeutidae", "Cossidae",
                                    "Saturniidae", "Choreutidae",
                                    "Scythrididae", "Lymantriidae",
                                    "Coleophoridae", "Notodontidae", "Lasiocampidae",
                                    "Eupterotidae","Elachistidae", "Euteliidae")
colnames(new.lepidopteraYesNo) <- c("Konye-wet", "Konye-dry", "Ayos-wet",
                                    "Ayos-dry", "Bokito-wet", "Bokito-dry")

library("pheatmap")
pheatmap(new.lepidopteraYesNo, color = c("indianred", "darkseagreen"), border_color = "white",
         cluster_rows = FALSE, cluster_cols = FALSE, breaks = c(0,0.5,1), angle_col = 0,
         legend = FALSE)
#####
