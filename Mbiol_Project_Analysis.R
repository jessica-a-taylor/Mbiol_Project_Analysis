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
filter.rows <- function(x, y, z) {                 
  z <- c()                                          
  for (row in 1:nrow(x)) {                         
    which.row <- which(x[row,] < y)                
    if (length(which.row)==length(x[row,])) {      
      z <- append(z, row) 
    }                                               
  } 
  return(z)
}
#####     
insect.data.cut <- filter.rows(insect.data[,2:270], 10, z)
                                                                   
insect.data <- insect.data[-c(insect.data.cut),]                
                                
# normalise to remove OTUs that are represented by less than 1% of reads in each sample
new.insect.data <- as.data.frame(insect.data[,-c(1,271:276)])
#####
   sum.reads <- c()                                                
   for (nid.sum in 1:ncol(new.insect.data)) {                      
     sum.reads <- append(sum.reads, sum(new.insect.data[,nid.sum]))
   }                                                               
                                                                   
   new.insect.data <- new.insect.data[,-c(247,248)]
   sum.reads <- c()                                                
   for (nid.sum in 1:ncol(new.insect.data)) {                      
     sum.reads <- append(sum.reads, sum(new.insect.data[,nid.sum]))
   }                                                               
                                                                   
   colnames(new.insect.data) <- c(sum.reads) 
                                                                   
   for (idc in as.character(sum.reads)) {                          
     normalised.diet <- (new.insect.data[,idc]/as.numeric(idc))*100
     new.insect.data[,idc] <- normalised.diet                      
   }                                                               
#####
   
new.insect.data.cut <- filter.rows(new.insect.data, 1, z)                                                                 
insect.data <- insect.data[-c(new.insect.data.cut),]            
insect.data <- insect.data[,-c(248,249)]

# diet calculations
new.insect.data <- as.data.frame(insect.data[,-c(1,269:274)])

# function for calculating FOO
#####
FOO.function <- function() {                                                                
                                                                                              insect.FOO.counts <- c()                                                                  
  for (seq_row in 1:nrow(new.insect.data)) {                                                
    I <- 0                                                                                  
    for (seq_col in 1:ncol(new.insect.data)) {                                              
      ifelse(new.insect.data[seq_row, seq_col] > 0,                                         
             I <- I + 1, I <- I) 
    }                                                                                       
    insect.FOO.counts <- append(insect.FOO.counts, (1/length(names(new.insect.data)))*I*100)
  }  
  return(insect.FOO.counts)
}
#####

insect.FOO.counts <- FOO.function()
insect.data <- cbind(insect.data, insect.FOO.counts)

# function for calculating POO
#####
POO.function <- function() {                                      
  insect.POO.counts <- c()                                        
  for (seq_row in 1:nrow(new.insect.data)) {                      
    I <- 0                                                        
    for (seq_col in 1:ncol(new.insect.data)) {                    
      ifelse(new.insect.data[seq_row, seq_col] > 0,               
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

insect.POO <- POO.function()
insect.data <- cbind(insect.data, insect.POO)

# function for calculating wPOO
#####
   wPOO.function <- function() {                                                                              
     for (row.no in 1:nrow(new.insect.data)) {                                                                
       for (col.no in 1:ncol(new.insect.data)) {                                                              
         if (new.insect.data[row.no, col.no] >0) {                                                            
           new.insect.data[row.no, col.no] <- 1
         }                                                                                                    
       }                                                                                                      
     }
                                                                                                              
     insect.sample.diet <- c()                                                                                
     for (seq_col in 1:ncol(new.insect.data)) {                                                               
     J <- 0
       for (seq_row in 1:nrow(new.insect.data)) {                                                             
         ifelse(new.insect.data[seq_row, seq_col] > 0, J <- J + 1, J <- J)                                    
       }                                                                                                      
       insect.sample.diet <- append(insect.sample.diet, J)                                                    
     }
     colnames(new.insect.data) <- c(insect.sample.diet)                                                       
     insect.value.df <- new.insect.data                                                                       
                                                                                                              
     new.insect.data <- as.data.frame(new.insect.data)                                                        
     reads <- c()
     for (insect.col in 1:ncol(new.insect.data)) {                                                            
       reads <- append(reads, new.insect.data[,insect.col])                                                   
     }                                                                                                        

     insect.df <- data.frame(Row = rep(1:944, times = 267),                                                   
                             Reads = reads,                                                                   
                             Total.reads = rep(insect.sample.diet, each = 944))                               

     insect.df <- cbind(insect.df, insect.df$Reads/insect.df$Total.reads)                                     
                                                                                                              
     insect.value.df <- new.insect.data                                                                       
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

insect.wPOO <- wPOO.function()
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

# flip columns and rows
zbj.focal <- t(zbj.focal)
zbj.focal <- as.data.frame(zbj.focal)
# function for converting seq data to occurrence data
##### 
seq.to.occurrence <- function(df) { 
  ifelse(df > 0,                    
  df <- 1, df <- 0)                 
}
#####

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
focal.info <- focal.info[,-c(1,2,5)]
colnames(zbj.focal) <- c(insect.data$OTU, "Lab.nbr", "Species", "Location", "Season")
#####

write_xlsx(zbj.focal, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\focal.species.data.xlsx")

# ordination plot for diet overlap between species
#####
library(vegan)                                                                              
library(dplyr)                                                                              
library(ggplot2)                                                                            
library(ggalt)                                                                              
library(ggforce)                                                                            
library(concaveman)                                                                         
                                                                                            
focal.nmds <- metaMDS(zbj.focal[1:944], distance = "jaccard")
plot(focal.nmds, type = "t", main = paste("NMDS/Jaccard - Stress = ",                       
                                             round(focal.nmds$stress, 3)))                     
                                                                                               
focal.envfit <- envfit(focal.nmds, focal.info, permutations = 999, na.rm = TRUE)            
plot(focal.envfit)                                                                          
                                                                                               
en_coord_cat <- as.data.frame(scores(focal.envfit, "factors")) 
                                                                                               
data.scores <- as.data.frame(scores(focal.nmds))
data.scores <- cbind(data.scores, focal.info)                                               
                                                                                               
data.scores <- data.scores[-c(155, 158),]

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

# db-rda for diet overlap between species
#####
focal.interact.species <- capscale(zbj.focal[1:944] ~ zbj.focal$Species, distance = "jaccard", add = TRUE)
anova(focal.interact.species, step = 1000, perm.max = 1000)

focal.interact.season <- capscale(zbj.focal[1:944] ~ zbj.focal$Season, distance = "jaccard", add = TRUE)
anova(focal.interact.season, step = 1000, perm.max = 1000)

focal.interact.location <- capscale(zbj.focal[1:944] ~ zbj.focal$Location, distance = "jaccard", add = TRUE)
anova(focal.interact.location, step = 1000, perm.max = 1000)

focal.interact.all <- capscale(zbj.focal[1:944] ~ Species + Season + Location, data = zbj.focal[946:948], distance = "jaccard", add = TRUE)
plot(focal.interact.all)
output <- anova(focal.interact.all, step = 1000, perm.max = 1000)

vec = output$SumOfSqs/sum(output$SumOfSqs)*100
table = output
table$SumOfSqs = vec
table

anova(focal.interact.all, by="axis", perm.max=500)

output2 <- anova(focal.interact.all, by="terms", permu=200)
vec2 = output2$SumOfSqs/sum(output2$SumOfSqs)*100
table2 = output2
table2$SumOfSqs = vec2
table2
#####

# ordination for diet overlap between seasons and locations for each species
#####
new.data.scores <- data.scores
new.data.scores <- cbind(new.data.scores, paste(new.data.scores$Location, "-", new.data.scores$Season))
colnames(new.data.scores) <- c("NMDS1", "NMDS2", "Species", "Location", "Season", "Location - Season")

new.data.scores$`Location - Season` <- factor(new.data.scores$`Location - Season`,
                                              levels = c("Konye - wet", "Konye - dry",
                                                         "Ayos - wet", "Ayos - dry",
                                                         "Bokito - wet", "Bokito - dry"))

colnames(new.data.scores) <- c("NMDS1", "NMDS2", "Species", "Location", "Season", "Location - Season")

LS.list <- unique(new.data.scores$`Location - Season`)
LS.list <- factor(LS.list, levels = c("Konye - wet", "Konye - dry",
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
#####

# db-rda for diet overlap between seasons and locations for each species
#####
new.zbj.focal <- zbj.focal
new.zbj.focal <- cbind(new.zbj.focal, paste(new.zbj.focal$Location, "-", new.zbj.focal$Season))
colnames(new.zbj.focal) <- c(names(zbj.focal), "Location - Season")

capscale.function <- function(df) {
  
  focal.interact <- capscale(df[1:944] ~ df$Species, distance = "jaccard", add = TRUE)
  interact.output <- anova(focal.interact, step = 1000, perm.max = 200)
  
  vec = interact.output$SumOfSqs/sum(interact.output$SumOfSqs)*100
  table = interact.output
  table$SumOfSqs = vec
  table

  return(table)
}

KonyeWet.rda <- capscale.function(new.zbj.focal[new.zbj.focal$`Location - Season`=="Konye - wet",]) 
KonyeDry.rda <- capscale.function(new.zbj.focal[new.zbj.focal$`Location - Season`=="Konye - dry",]) 
AyosWet.rda <- capscale.function(new.zbj.focal[new.zbj.focal$`Location - Season`=="Ayos - wet",]) 
AyosDry.rda <- capscale.function(new.zbj.focal[new.zbj.focal$`Location - Season`=="Ayos - dry",]) 
BokitoWet.rda <- capscale.function(new.zbj.focal[new.zbj.focal$`Location - Season`=="Bokito - wet",]) 
BokitoDry.rda <- capscale.function(new.zbj.focal[new.zbj.focal$`Location - Season`=="Bokito - dry",]) 
#####

# total number of reads for all species == OTU.sum
#####
insect.data.reads <- insect.data[,-c(1,269:277)]          
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
                                                 
RA.nbrs <- zbj.insect.data.names[RA.focal.nbrs]
RA.focal.df <- zbj.insect.data[,c(RA.nbrs)]
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

# filter out OTUs with 0 occurrences in the focal datasets 
#####
HR.focal.df.cut <- filter.rows(HR.focal.df[1:79], 1, z)
HR.focal.df <- HR.focal.df[-c(HR.focal.df.cut),]

RA.focal.df.cut <- filter.rows(RA.focal.df[1:82], 1, z)
RA.focal.df <- RA.focal.df[-c(RA.focal.df.cut),]
#####

# taxonomic resolution
#####
insect.data.order <- insect.data[,c(1,271:274,277)]

HR.focal.df <- cbind(HR.focal.df, insect.data.order[-c(HR.focal.df.cut),])
RA.focal.df <- cbind(RA.focal.df, insect.data.order[-c(RA.focal.df.cut),])

#   not.na <- c()                            
#   for (row in 1:nrow(RA.focal.df)) {       
#     ifelse(RA.focal.df[row, "genus"]=="NA",
#            not.na <- append(not.na, row),  
#            not.na <- not.na)               
#   }                                        
#                                            
#   length(RA.focal.df$order)-length(not.na) 
#   (90/length(RA.focal.df$order))*100
#####

# convert to occurrence data
#####
HR.focal.df <- data.frame(lapply(HR.focal.df[,1:79],seq.to.occurrence))
RA.focal.df <- data.frame(lapply(RA.focal.df[,1:82],seq.to.occurrence))

HR.focal.df <- cbind(HR.focal.df, insect.data.order[-c(HR.focal.df.cut),])
RA.focal.df <- cbind(RA.focal.df, insect.data.order[-c(RA.focal.df.cut),])

colnames(HR.focal.df) <- c(HR.nbrs, "OTU", "order", "family", "genus", "species", "insect.wPOO")
colnames(RA.focal.df) <- c(RA.nbrs, "OTU", "order", "family", "genus", "species", "insect.wPOO")
#####

# Venn diagrams showing exclusive and shared OTUs consumed by the focal species
#####
library(RAM)
focal.compare.OTUs <- list(H.ruber = HR.focal.df$OTU,
                           R.alcyone = RA.focal.df$OTU)

group.venn(focal.compare.OTUs, label=FALSE,  lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(330, 150), cat.dist = 0.05)
#####

# Venn diagrams showing exclusive and shared OTUs consumed by the focal species 
#     in each location
#####
HR.zbj.data <- zbj_final[zbj_final$Species=="Hipposideros ruber",]
RA.zbj.data <- zbj_final[zbj_final$Species=="Rhinolophus alcyone",]

HR.Konye.data <- HR.zbj.data[HR.zbj.data$Location=="Konye",]
HR.Konye.data <- HR.Konye.data[-c(which(is.na(HR.Konye.data$Lab.nbr))),]
HR.Konye.nbrs <- HR.Konye.data$Lab.nbr

HR.Konye.focal.df <- HR.focal.df[,c(as.character(HR.Konye.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.Konye.focal.df.cut <- filter.rows(HR.Konye.focal.df[,1:37], 1, z)
HR.Konye.focal.df <- HR.Konye.focal.df[-c(HR.Konye.focal.df.cut),]

HR.Ayos.data <- HR.zbj.data[HR.zbj.data$Location=="Ayos",]
HR.Ayos.data <- HR.Ayos.data[-c(which(is.na(HR.Ayos.data$Lab.nbr))),]
HR.Ayos.nbrs <- HR.Ayos.data$Lab.nbr

HR.Ayos.focal.df <- HR.focal.df[,c(as.character(HR.Ayos.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.Ayos.focal.df.cut <- filter.rows(HR.Ayos.focal.df[,1:35], 1, z)
HR.Ayos.focal.df <- HR.Ayos.focal.df[-c(HR.Ayos.focal.df.cut),]

HR.Bokito.data <- HR.zbj.data[HR.zbj.data$Location=="Bokito",]
HR.Bokito.data <- HR.Bokito.data[-c(which(is.na(HR.Bokito.data$Lab.nbr))),]
HR.Bokito.nbrs <- HR.Bokito.data$Lab.nbr

HR.Bokito.focal.df <- HR.focal.df[,c(as.character(HR.Bokito.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.Bokito.focal.df.cut <- filter.rows(HR.Bokito.focal.df[,1:7], 1, z)
HR.Bokito.focal.df <- HR.Bokito.focal.df[-c(HR.Bokito.focal.df.cut),]

RA.Konye.data <- RA.zbj.data[RA.zbj.data$Location=="Konye",]
RA.Konye.data <- RA.Konye.data[-c(which(is.na(RA.Konye.data$Lab.nbr))),]
RA.Konye.nbrs <- RA.Konye.data$Lab.nbr

RA.Konye.focal.df <- RA.focal.df[,c(as.character(RA.Konye.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.Konye.focal.df.cut <- filter.rows(RA.Konye.focal.df[,1:47], 1, z)
RA.Konye.focal.df <- RA.Konye.focal.df[-c(RA.Konye.focal.df.cut),]

RA.Ayos.data <- RA.zbj.data[RA.zbj.data$Location=="Ayos",]
RA.Ayos.data <- RA.Ayos.data[-c(which(is.na(RA.Ayos.data$Lab.nbr))),]
RA.Ayos.nbrs <- RA.Ayos.data$Lab.nbr

RA.Ayos.focal.df <- RA.focal.df[,c(as.character(RA.Ayos.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.Ayos.focal.df.cut <- filter.rows(RA.Ayos.focal.df[,1:33], 1, z)
RA.Ayos.focal.df <- RA.Ayos.focal.df[-c(RA.Ayos.focal.df.cut),]

RA.Bokito.data <- RA.zbj.data[RA.zbj.data$Location=="Bokito",]
RA.Bokito.data <- RA.Bokito.data[-c(which(is.na(RA.Bokito.data$Lab.nbr))),]
RA.Bokito.nbrs <- RA.Bokito.data$Lab.nbr

RA.Bokito.focal.df <- RA.focal.df[,c(as.character(RA.Bokito.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.Bokito.focal.df.cut <- filter.rows(RA.Bokito.focal.df[,1:2], 1, z)
RA.Bokito.focal.df <- RA.Bokito.focal.df[-c(RA.Bokito.focal.df.cut),]

focal.compare.locations.OTUs <- list(H.ruber = HR.Konye.focal.df$OTU,
                                     R.alcyone = RA.Konye.focal.df$OTU)

group.venn(focal.compare.locations.OTUs, label=FALSE,  lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(400, 200), cat.dist = 0.05)
#####

# Venn diagrams showing exclusive and shared OTUs consumed by the focal species 
#     in each season
#####
HR.wet.data <- HR.zbj.data[HR.zbj.data$Season=="wet",]
HR.wet.data <- HR.wet.data[-c(which(is.na(HR.wet.data$Lab.nbr))),]
HR.wet.nbrs <- HR.wet.data$Lab.nbr

HR.wet.focal.df <- HR.focal.df[,c(as.character(HR.wet.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.wet.focal.df.cut <- filter.rows(HR.wet.focal.df[,1:7], 1, z)
HR.wet.focal.df <- HR.wet.focal.df[-c(HR.wet.focal.df.cut),]

RA.wet.data <- RA.zbj.data[RA.zbj.data$Season=="wet",]
RA.wet.data <- RA.wet.data[-c(which(is.na(RA.wet.data$Lab.nbr))),]
RA.wet.nbrs <- RA.wet.data$Lab.nbr

RA.wet.focal.df <- RA.focal.df[,c(as.character(RA.wet.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.wet.focal.df.cut <- filter.rows(RA.wet.focal.df[,1:7], 1, z)
RA.wet.focal.df <- RA.wet.focal.df[-c(RA.wet.focal.df.cut),]

HR.dry.data <- HR.zbj.data[HR.zbj.data$Season=="dry",]
HR.dry.data <- HR.dry.data[-c(which(is.na(HR.dry.data$Lab.nbr))),]
HR.dry.nbrs <- HR.dry.data$Lab.nbr

HR.dry.focal.df <- HR.focal.df[,c(as.character(HR.dry.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.dry.focal.df.cut <- filter.rows(HR.dry.focal.df[,1:7], 1, z)
HR.dry.focal.df <- HR.dry.focal.df[-c(HR.dry.focal.df.cut),]

RA.dry.data <- RA.zbj.data[RA.zbj.data$Season=="dry",]
RA.dry.data <- RA.dry.data[-c(which(is.na(RA.dry.data$Lab.nbr))),]
RA.dry.nbrs <- RA.dry.data$Lab.nbr

RA.dry.focal.df <- RA.focal.df[,c(as.character(RA.dry.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.dry.focal.df.cut <- filter.rows(RA.dry.focal.df[,1:7], 1, z)
RA.dry.focal.df <- RA.dry.focal.df[-c(RA.dry.focal.df.cut),]

focal.compare.seasons.OTUs <- list(H.ruber = HR.wet.focal.df$OTU,
                                   R.alcyone = RA.wet.focal.df$OTU)

group.venn(focal.compare.seasons.OTUs, label=FALSE,  lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(330, 150), cat.dist = 0.05)
#####

# Venn diagrams showing exclusive and shared OTUs consumed by the focal species
#     in each location during each season
#####
new.HR.zbj.data <- HR.zbj.data
new.HR.zbj.data <- cbind(new.HR.zbj.data, paste(new.HR.zbj.data$Location, "-" ,new.HR.zbj.data$Season))
colnames(new.HR.zbj.data) <- c(names(HR.zbj.data), "Location - Season")

HR.KW.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Konye - wet",]
HR.KW.nbrs <- HR.KW.data$Lab.nbr

HR.KW.focal.df <- HR.focal.df[,c(as.character(HR.KW.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.KW.focal.df.cut <- filter.rows(HR.KW.focal.df[,1:13], 1, z)
HR.KW.focal.df <- HR.KW.focal.df[-c(HR.KW.focal.df.cut),]

HR.KD.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Konye - dry",]
HR.KD.nbrs <- HR.KD.data$Lab.nbr

HR.KD.focal.df <- HR.focal.df[,c(as.character(HR.KD.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.KD.focal.df.cut <- filter.rows(HR.KD.focal.df[,1:24], 1, z)
HR.KD.focal.df <- HR.KD.focal.df[-c(HR.KD.focal.df.cut),]

HR.AW.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Ayos - wet",]
HR.AW.nbrs <- HR.AW.data$Lab.nbr

HR.AW.focal.df <- HR.focal.df[,c(as.character(HR.AW.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.AW.focal.df.cut <- filter.rows(HR.AW.focal.df[,1:10], 1, z)
HR.AW.focal.df <- HR.AW.focal.df[-c(HR.AW.focal.df.cut),]

HR.AD.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Ayos - dry",]
HR.AD.nbrs <- HR.AD.data$Lab.nbr

HR.AD.focal.df <- HR.focal.df[,c(as.character(HR.AD.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.AD.focal.df.cut <- filter.rows(HR.AD.focal.df[,1:25], 1, z)
HR.AD.focal.df <- HR.AD.focal.df[-c(HR.AD.focal.df.cut),]

HR.BW.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Bokito - wet",]
HR.BW.nbrs <- HR.BW.data$Lab.nbr

HR.BW.focal.df <- HR.focal.df[,c(as.character(HR.BW.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.BW.focal.df.cut <- filter.rows(HR.BW.focal.df[,1:4], 1, z)
HR.BW.focal.df <- HR.BW.focal.df[-c(HR.BW.focal.df.cut),]

HR.BD.data <- new.HR.zbj.data[new.HR.zbj.data$`Location - Season`=="Bokito - dry",]
HR.BD.nbrs <- HR.BD.data$Lab.nbr

HR.BD.focal.df <- HR.focal.df[,c(as.character(HR.BD.nbrs), "OTU", "order", "family", "insect.wPOO")]
HR.BD.focal.df.cut <- filter.rows(HR.BD.focal.df[,1:3], 1, z)
HR.BD.focal.df <- HR.BD.focal.df[-c(HR.BD.focal.df.cut),]

new.RA.zbj.data <- RA.zbj.data
new.RA.zbj.data <- cbind(new.RA.zbj.data, paste(new.RA.zbj.data$Location, "-" ,new.RA.zbj.data$Season))
colnames(new.RA.zbj.data) <- c(names(RA.zbj.data), "Location - Season")

RA.KW.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Konye - wet",]
RA.KW.nbrs <- RA.KW.data$Lab.nbr

RA.KW.focal.df <- RA.focal.df[,c(as.character(RA.KW.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.KW.focal.df.cut <- filter.rows(RA.KW.focal.df[,1:29], 1, z)
RA.KW.focal.df <- RA.KW.focal.df[-c(RA.KW.focal.df.cut),]

RA.KD.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Konye - dry",]
RA.KD.nbrs <- RA.KD.data$Lab.nbr

RA.KD.focal.df <- RA.focal.df[,c(as.character(RA.KD.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.KD.focal.df.cut <- filter.rows(RA.KD.focal.df[,1:18], 1, z)
RA.KD.focal.df <- RA.KD.focal.df[-c(RA.KD.focal.df.cut),]

RA.AW.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Ayos - wet",]
RA.AW.nbrs <- RA.AW.data$Lab.nbr

RA.AW.focal.df <- RA.focal.df[,c(as.character(RA.AW.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.AW.focal.df.cut <- filter.rows(RA.AW.focal.df[,1:20], 1, z)
RA.AW.focal.df <- RA.AW.focal.df[-c(RA.AW.focal.df.cut),]

RA.AD.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Ayos - dry",]
RA.AD.nbrs <- RA.AD.data$Lab.nbr

RA.AD.focal.df <- RA.focal.df[,c(as.character(RA.AD.nbrs), "OTU", "order", "family", "insect.wPOO")]
RA.AD.focal.df.cut <- filter.rows(RA.AD.focal.df[,1:13], 1, z)
RA.AD.focal.df <- RA.AD.focal.df[-c(RA.AD.focal.df.cut),]

RA.BW.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Bokito - wet",]
RA.BW.nbrs <- RA.BW.data$Lab.nbr

RA.BW.focal.df <- RA.focal.df[,c(as.character(RA.BW.nbrs), "OTU", "order", "family", "insect.wPOO")]

RA.BD.data <- new.RA.zbj.data[new.RA.zbj.data$`Location - Season`=="Bokito - dry",]
RA.BD.nbrs <- RA.BD.data$Lab.nbr

RA.BD.focal.df <- RA.focal.df[,c(as.character(RA.BD.nbrs), "OTU", "order", "family", "insect.wPOO")]


focal.compare.LS.OTUs <- list(H.ruber = HR.BD.focal.df$OTU,
                                   R.alcyone = RA.BD.focal.df$OTU)

group.venn(focal.compare.LS.OTUs, label=FALSE,  lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.pos=c(330, 150), cat.dist = 0.05)
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

# compare average number of OTUs between focal species in each location
#####
HR.count.OTUs <- count.OTUs(HR.wet.focal.df[,1:27])
RA.count.OTUs <- count.OTUs(RA.wet.focal.df[,1:50])

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

# bootstrapping function to test for differences between medians (non-parametric test)
#####
library(boot)
med.diff <- function(a,b) {                         
  df <- a[b,]                                    
  median(df$Count[df$Species=="H.ruber"]) -  
    median(df$Count[df$Species=="R.alcyone"])
}

boot.out <- boot(data = focal.compare.OTUs.df, statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
#####

# filter for OTUs assigned to taxonimic order
#####
HR.na <- which(HR.focal.df$order=="NA", arr.ind = TRUE)      
HR.focal.df <- HR.focal.df[-c(HR.na),]

RA.na <- which(RA.focal.df$order=="NA", arr.ind = TRUE)      
RA.focal.df <- RA.focal.df[-c(RA.na),]

(length(which(HR.focal.df$order=="Lepidoptera"))/length(HR.focal.df$order))*100
(length(which(RA.focal.df$order=="Lepidoptera"))/length(RA.focal.df$order))*100

#####

# filter for OTUs assigned to taxonomic family
#####
HR.na <- which(HR.focal.df$family=="NA", arr.ind = TRUE)      
HR.focal.df <- HR.focal.df[-c(HR.na),]
                                                                               
RA.na <- which(RA.focal.df$family=="NA", arr.ind = TRUE)      
RA.focal.df <- RA.focal.df[-c(RA.na),]
#####

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
#####

# diet composition - most common Lepidopteran families in each location and season
#####
HR.Lepidoptera <- HR.focal.df[HR.focal.df$order=="Lepidoptera",]
HR.Lepidoptera.fam <- length(unique(HR.Lepidoptera$family))

RA.wet.Lepidoptera <- RA.wet.focal.df[RA.wet.focal.df$order=="Lepidoptera",]
RA.wet.Lepidoptera.fam <- unique(RA.wet.Lepidoptera$family)

RA.dry.Lepidoptera <- RA.dry.focal.df[RA.dry.focal.df$order=="Lepidoptera",]
RA.dry.Lepidoptera.fam <- unique(RA.dry.Lepidoptera$family)

RA.Bokito.Lepidoptera <- RA.Bokito.focal.df[RA.Bokito.focal.df$order=="Lepidoptera",]
RA.Bokito.Lepidoptera.fam <- unique(RA.Bokito.Lepidoptera$family)

RA.Lepidopteran.location.fam <- list(Wet = RA.wet.Lepidoptera.fam,
                                  Dry = RA.dry.Lepidoptera.fam)

group.venn(RA.Lepidopteran.location.fam, label=FALSE,  lab.cex=1, cex = 1.5,
           cat.cex = 1.5, cat.dist = 0.05, cat.pos = c(320, 400))
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

HR.Shannon.index <- Shannon.index.function(HR.focal.df[,1:79])
RA.Shannon.index <- Shannon.index.function(RA.focal.df[,1:82])

shapiro.test(c(HR.Shannon.index, RA.Shannon.index))

focal.compare.Shannon <- data.frame(Species = c(rep("H.ruber", times = length(HR.Shannon.index)),
                                                rep("R.alcyone", times = length(RA.Shannon.index))),
                                    Count = c(HR.Shannon.index, RA.Shannon.index))

focal.compare.Shannon.boxplot <- ggplot(focal.compare.Shannon, aes(x = Species, y = Count, fill = Species)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                axis.title.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 12),
                                                axis.title.y = element_text(size = 14),
                                                axis.text.x = element_text(size = 12)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Species", y = "Shannon index") 


wilcox.test(HR.Shannon.index, RA.Shannon.index)
boot.out <- boot(data = focal.compare.Shannon, statistic = med.diff, R = 1000)
median(boot.out$t)

boot.ci(boot.out, type = "perc")

# Shannon index for each location
#####
HR.Konye.Shannon.index <- Shannon.index.function(HR.Konye.focal.df[,1:37])
HR.Ayos.Shannon.index <- Shannon.index.function(HR.Ayos.focal.df[,1:35])
HR.Bokito.Shannon.index <- Shannon.index.function(HR.Bokito.focal.df[,1:7])

RA.Konye.Shannon.index <- Shannon.index.function(RA.Konye.focal.df[,1:47])
RA.Ayos.Shannon.index <- Shannon.index.function(RA.Ayos.focal.df[,1:33])
RA.Bokito.Shannon.index <- Shannon.index.function(RA.Bokito.focal.df[,1:2])

locations.compare.Shannon <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.Konye.Shannon.index, HR.Ayos.Shannon.index, HR.Bokito.Shannon.index))),
                                                       rep("R.alcyone", times = length(c(RA.Konye.Shannon.index, RA.Ayos.Shannon.index, RA.Bokito.Shannon.index)))),
                                           Location = c(rep("Konye", times = length(HR.Konye.Shannon.index)),
                                                        rep("Ayos", times = length(HR.Ayos.Shannon.index)),
                                                        rep("Bokito", times = length(HR.Bokito.Shannon.index)),
                                                        rep("Konye", times = length(RA.Konye.Shannon.index)),
                                                        rep("Ayos", times = length(RA.Ayos.Shannon.index)),
                                                        rep("Bokito", times = length(RA.Bokito.Shannon.index))),
                                           Shannon.index = c(HR.Konye.Shannon.index, HR.Ayos.Shannon.index, HR.Bokito.Shannon.index,
                                                             RA.Konye.Shannon.index, RA.Ayos.Shannon.index, RA.Bokito.Shannon.index))

location.compare.Shannon.boxplot <- ggplot(locations.compare.Shannon, aes(x = Location, y = Shannon.index, fill = Location)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 14),
                                                axis.title.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 12),
                                                axis.title.y = element_text(size = 14),
                                                axis.text.x = element_text(size = 12)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Species", y = "Shannon index") + facet_wrap(vars(Species))   


HR.Shannon <- locations.compare.Shannon[locations.compare.Shannon$Species=="H.ruber",]
RA.Shannon <- locations.compare.Shannon[locations.compare.Shannon$Species=="R.alcyone",]

wilcox.test(HR.Shannon$Shannon.index, RA.Shannon$Shannon.index)

Shannon.kruskal <- kruskal.test(data = locations.compare.Shannon, x = locations.compare.Shannon$Shannon.index,
             g = locations.compare.Shannon$Location, formula = Shannon.index ~ Location)



med.diff <- function(d, i) {                         
  df.shannon <- d[i,]                                    
  median(df.shannon$Shannon.index[df.shannon$Species=="H.ruber"]) -  
    median(df.shannon$Shannon.index[df.shannon$Species=="R.alcyone"])
}

boot.out <- boot(data = locations.compare.Shannon[,-2], 
                 statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
#####

# Shannon index for each seasons
#####
HR.wet.Shannon.index <- Shannon.index.function(HR.wet.focal.df[,1:27])
HR.dry.Shannon.index <- Shannon.index.function(HR.dry.focal.df[,1:52])

RA.wet.Shannon.index <- Shannon.index.function(RA.wet.focal.df[,1:50])
RA.dry.Shannon.index <- Shannon.index.function(RA.dry.focal.df[,1:32])

seasons.compare.Shannon <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.wet.Shannon.index, HR.dry.Shannon.index))),
                                                    rep("R.alcyone", times = length(c(RA.wet.Shannon.index, RA.dry.Shannon.index)))),
                                        Season = c(rep("Wet", times = length(HR.wet.Shannon.index)),
                                                     rep("Dry", times = length(HR.dry.Shannon.index)),
                                                     rep("Wet", times = length(RA.wet.Shannon.index)),
                                                     rep("Dry", times = length(RA.dry.Shannon.index))),
                                        Shannon.index = c(HR.wet.Shannon.index, HR.dry.Shannon.index, 
                                                          RA.wet.Shannon.index, RA.dry.Shannon.index))

seasons.compare.Shannon.boxplot <- ggplot(seasons.compare.Shannon, aes(x = Season, y = Shannon.index, fill = Season)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 14),
                                                axis.title.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 12),
                                                axis.title.y = element_text(size = 14),
                                                axis.text.x = element_text(size = 12)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Species", y = "Shannon index") + facet_wrap(vars(Species))   


wet.Shannon <- seasons.compare.Shannon[seasons.compare.Shannon$Season=="Wet",]
dry.Shannon <- seasons.compare.Shannon[seasons.compare.Shannon$Season=="Dry",]

Shannon.wilcox <- wilcox.test(wet.Shannon$Shannon.index, dry.Shannon$Shannon.index)

boot.out <- boot(data = seasons.compare.Shannon[,-2], 
                 statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
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
HR.Levins.index <- Levins.function(HR.focal.df[,1:79])
RA.Levins.index <- Levins.function(RA.focal.df[,1:82])

shapiro.test(c(HR.Levins.index, RA.Levins.index))

wilcox.test(HR.Shannon.index, RA.Shannon.index)
boot.out <- boot(data = focal.compare.Shannon, statistic = med.diff, R = 1000)
median(boot.out$t)

boot.ci(boot.out, type = "perc")

HR.Konye.Levins <- Levins.function(HR.Konye.focal.df[,1:37])
HR.Ayos.Levins <- Levins.function(HR.Ayos.focal.df[,1:35])
HR.Bokito.Levins <- Levins.function(HR.Bokito.focal.df[,1:7])

RA.Konye.Levins <- Levins.function(RA.Konye.focal.df[,1:47])
RA.Ayos.Levins <- Levins.function(RA.Ayos.focal.df[,1:33])
RA.Bokito.Levins <- Levins.function(RA.Bokito.focal.df[,1:2])

HR.wet.Levins <- Levins.function(HR.wet.focal.df[,1:27])
HR.dry.Levins <- Levins.function(HR.dry.focal.df[,1:52])

RA.wet.focal.df <- RA.wet.focal.df[,-20] # no occurrences in this column
RA.dry.focal.df <- RA.dry.focal.df[,-13] # no occurrences in this column
RA.wet.Levins <- Levins.function(RA.wet.focal.df[,1:49])
RA.dry.Levins <- Levins.function(RA.dry.focal.df[,1:31])


# Levin's index for each location
#####
locations.compare.Levin <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.Konye.Levins$All.Levins, HR.Ayos.Levins$All.Levins, HR.Bokito.Levins$All.Levins))),
                                                  rep("R.alcyone", times = length(c(RA.Konye.Levins$All.Levins, RA.Ayos.Levins$All.Levins, RA.Bokito.Levins$All.Levins)))),
                                      Location = c(rep("Konye", times = length(HR.Konye.Levins$All.Levins)),
                                                 rep("Ayos", times = length(HR.Ayos.Levins$All.Levins)),
                                                 rep("Bokito", times = length(HR.Bokito.Levins$All.Levins)),
                                                 rep("Konye", times = length(RA.Konye.Levins$All.Levins)),
                                                 rep("Ayos", times = length(RA.Ayos.Levins$All.Levins)),
                                                 rep("Bokito", times = length(RA.Bokito.Levins$All.Levins))),
                                      Levins = c(HR.Konye.Levins$All.Levins.std, HR.Ayos.Levins$All.Levins.std, HR.Bokito.Levins$All.Levins.std, 
                                                      RA.Konye.Levins$All.Levins.std, RA.Ayos.Levins$All.Levins.std, RA.Bokito.Levins$All.Levins.std))

y_expression <- expression("Levin's index" ~ "" ~(B[A]))

library(ggpubr)
locations.compare.Levin.boxplot <- ggplot(locations.compare.Levin, aes(x = Location, y = Levins, fill = Location)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 14),
                                                axis.title.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 12),
                                                axis.title.y = element_text(size = 14),
                                                axis.text.x = element_text(size = 12)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location", y = y_expression) + facet_wrap(vars(Species)) 

shapiro.test(locations.compare.Levin$Levins)
HR.Levins <- locations.compare.Levin[locations.compare.Levin$Species=="H.ruber",]
RA.Levins <- locations.compare.Levin[locations.compare.Levin$Species=="R.alcyone",]

Levins.kruskal <- kruskal.test(data = locations.compare.Levin, x = locations.compare.Levin$Levins,
                                g = locations.compare.Levin$Location, formula = Levins ~ Location)

library(FSA)
dunnTest(Levins ~ Location, data = locations.compare.Levin, method = "bonferroni")

location.Levins <- locations.compare.Levin[locations.compare.Levin$Location=="Bokito",]

wilcox.test(HR.Levins$Levins, RA.Levins$Levins)

med.diff <- function(d, i) {                         
  df.levins <- d[i,]                                    
  median(df.levins$Levins[df.levins$Species=="H.ruber"]) -  
    median(df.levins$Levins[df.levins$Species=="R.alcyone"])
}

boot.out <- boot(data = location.Levins[,-2], 
                 statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
#####

# Levin's index for each season
#####
seasons.compare.Levin <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.wet.Levins$All.Levins, HR.dry.Levins$All.Levins))),
                                                  rep("R.alcyone", times = length(c(RA.wet.Levins$All.Levins, RA.dry.Levins$All.Levins)))),
                                      Season = c(rep("Wet", times = length(HR.wet.Levins$All.Levins)),
                                                   rep("Dry", times = length(HR.dry.Levins$All.Levins)),
                                                   rep("Wet", times = length(RA.wet.Levins$All.Levins)),
                                                   rep("Dry", times = length(RA.dry.Levins$All.Levins))),
                                      Levins = c(HR.wet.Levins$All.Levins.std, HR.dry.Levins$All.Levins.std, 
                                                 RA.wet.Levins$All.Levins.std, RA.dry.Levins$All.Levins.std))

seasons.compare.Levin.boxplot <- ggplot(seasons.compare.Levin, aes(x = Season, y = Levins, fill = Season)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 14),
                                                axis.title.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 12),
                                                axis.title.y = element_text(size = 14),
                                                axis.text.x = element_text(size = 12)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Season", y = y_expression) + facet_wrap(vars(Species)) 

shapiro.test(seasons.compare.Levin$Levins)

wet.Levins <- seasons.compare.Levin[seasons.compare.Levin$Season=="Wet",]
dry.Levins <- seasons.compare.Levin[seasons.compare.Levin$Season=="Dry",]

Levins.wilcox <- wilcox.test(wet.Levins$Levins, dry.Levins$Levins)

boot.out <- boot(data = dry.Levins[,-2], 
                 statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
#####

# Niche overlap
# function for calculating Pianka's index
#####
Pianka.function <- function(df) {
  HR.diet.average <- c()
  for (row in 1:nrow(df)) {
    all.diets <- c()
    for (col in 1:ncol(df)) {
      all.diets <- append(all.diets, df[row,col])
    }
    HR.diet.average <- append(HR.diet.average, mean(all.diets))
  }
  
  HR.diet.prop.Konye <- c()
  for (prop in HR.diet.average) {
    HR.diet.prop.Konye <- append(HR.diet.prop.Konye, prop/sum(HR.diet.average))
  }
  return(HR.diet.prop.Konye)
}
#####

HR.Konye.Pianka <- Pianka.function(HR.Konye.focal.df[,1:37])
HR.Ayos.Pianka <- Pianka.function(HR.Ayos.focal.df[,1:35])
HR.Bokito.Pianka <- Pianka.function(HR.Bokito.focal.df[,1:7])

RA.Konye.Pianka <- Pianka.function(RA.Konye.focal.df[,1:47])
RA.Ayos.Pianka <- Pianka.function(RA.Ayos.focal.df[,1:33])
RA.Bokito.Pianka <- Pianka.function(RA.Bokito.focal.df[,1:2])

HR.wet.Pianka <- Pianka.function(HR.wet.focal.df[,1:27])
HR.dry.Pianka <- Pianka.function(HR.dry.focal.df[,1:52])

RA.wet.Pianka <- Pianka.function(RA.wet.focal.df[,1:49])
RA.dry.Pianka <- Pianka.function(RA.dry.focal.df[,1:31])

# Pianka's index for each location
#####
locations.compare.Pianka <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.Konye.Pianka, HR.Ayos.Pianka, HR.Bokito.Pianka))),
                                                  rep("R.alcyone", times = length(c(RA.Konye.Pianka, RA.Ayos.Pianka, RA.Bokito.Pianka)))),
                                      Location = c(rep("Konye", times = length(HR.Konye.Pianka)),
                                                   rep("Ayos", times = length(HR.Ayos.Pianka)),
                                                   rep("Bokito", times = length(HR.Bokito.Pianka)),
                                                   rep("Konye", times = length(RA.Konye.Pianka)),
                                                   rep("Ayos", times = length(RA.Ayos.Pianka)),
                                                   rep("Bokito", times = length(RA.Bokito.Pianka))),
                                      Pianka = c(HR.Konye.Pianka, HR.Ayos.Pianka, HR.Bokito.Pianka, 
                                                 RA.Konye.Pianka, RA.Ayos.Pianka, RA.Bokito.Pianka))

y_expression <- expression("Pianka's index" ~ "" ~(O[jk]))

locations.compare.Pianka.boxplot <- ggplot(locations.compare.Pianka, aes(x = Location, y = Pianka, fill = Location)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Location", y = y_expression) + facet_wrap(vars(Species))

shapiro.test(locations.compare.Pianka$Pianka)

HR.Pianka <- locations.compare.Pianka[locations.compare.Pianka$Species=="H.ruber",]
RA.Pianka <- locations.compare.Pianka[locations.compare.Pianka$Species=="R.alcyone",]

wilcox.test(HR.Pianka$Pianka, RA.Pianka$Pianka)

Pianka.kruskal <- kruskal.test(data = locations.compare.Pianka, x = locations.compare.Pianka$Pianka,
                               g = locations.compare.Pianka$Location, formula = Pianka ~ Location)

dunnTest(Pianka ~ Location, data = locations.compare.Pianka, method = "bonferroni")

location.Pianka <- locations.compare.Pianka[locations.compare.Pianka$Location=="Ayos",]

med.diff <- function(d, i) {                         
  df.Pianka <- d[i,]                                    
  median(df.Pianka$Pianka[df.Pianka$Species=="H.ruber"]) -  
    median(df.Pianka$Pianka[df.Pianka$Species=="R.alcyone"])
}

boot.out <- boot(data = location.Pianka[,-2], 
                 statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
#####

# Pianka's index for each season
#####
seasons.compare.Pianka <- data.frame(Species = c(rep("H.ruber", times = length(c(HR.wet.Pianka, HR.dry.Pianka))),
                                                   rep("R.alcyone", times = length(c(RA.wet.Pianka, RA.dry.Pianka)))),
                                       Season = c(rep("Wet", times = length(HR.wet.Pianka)),
                                                    rep("Dry", times = length(HR.dry.Pianka)),
                                                    rep("Wet", times = length(RA.wet.Pianka)),
                                                    rep("Dry", times = length(RA.dry.Pianka))),
                                       Pianka = c(HR.wet.Pianka, HR.dry.Pianka, 
                                                  RA.wet.Pianka, RA.dry.Pianka))

y_expression <- expression("Pianka's index" ~ "" ~(O[jk]))

seasons.compare.Pianka.boxplot <- ggplot(seasons.compare.Pianka, aes(x = Season, y = Pianka, fill = Season)) +
  geom_boxplot(show.legend = FALSE) + theme_classic() +
  scale_fill_brewer(palette = "Accent") + theme(strip.background = element_blank(),
                                                strip.text = element_text(size = 12),
                                                axis.title.x = element_text(size = 12),
                                                axis.text.y = element_text(size = 10),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10)) +
  stat_summary(fun="mean", show.legend = FALSE, color="black", shape=18) + 
  labs(x = "Season", y = y_expression) + facet_wrap(vars(Species))

wet.Pianka <- seasons.compare.Pianka[seasons.compare.Pianka$Season=="Wet",]
dry.Pianka <- seasons.compare.Pianka[seasons.compare.Pianka$Season=="Dry",]

Pianka.wilcox <- wilcox.test(wet.Pianka$Pianka, dry.Pianka$Pianka)

boot.out <- boot(data = dry.Pianka[,-2], 
                 statistic = med.diff, R = 500)
median(boot.out$t)

boot.ci(boot.out, type = "perc")
#####

# Occurrence of pest sequences
pest.data <- readxl::read_xlsx("C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\Pest sequences.xlsx")
pest.data <- as.data.frame(pest.data)
pest.data <- pest.data[,-c(56,185,272,274:291)]


rownames(pest.data) <- c(pest.data$OTU)
rownames(HR.focal.df) <- c(HR.focal.df$OTU)
rownames(RA.focal.df) <- c(RA.focal.df$OTU)

focal.OTUs <- unique(c(HR.focal.df$OTU, RA.focal.df$OTU))
focal.OTUs <- focal.OTUs[which(focal.OTUs %in% pest.data$OTU)]
pest.OTUs <- pest.data[c(6,13,26,27,57,62,73,80:82,91,97:102,
                         111:114,121,129:137,142,149,170:183,
                         197:202,207:209,15,74,75,88,107,122,
                         150,203,210,211), "OTU"]

focal.pest.data <- pest.data[c(focal.OTUs),]

final.pest.data <- seq.to.occurrence(focal.pest.data[,2:270])
final.pest.data <- cbind(final.pest.data, focal.pest.data[,c(1, 273:276)])

HR.pest.df <- final.pest.data[,c(HR.nbrs)]
HR.pest.df <- cbind(HR.pest.df, final.pest.data[,c(270:274)])
HR.pest.df.cut <- filter.rows(HR.pest.df[,1:79], 1, z)
HR.pest.df <- HR.pest.df[-c(HR.pest.df.cut),]

RA.pest.df <- final.pest.data[,c(RA.nbrs)]
RA.pest.df <- cbind(RA.pest.df, final.pest.data[,c(270:274)])
RA.pest.df.cut <- filter.rows(RA.pest.df[1:82], 1, z)
RA.pest.df <- RA.pest.df[-c(RA.pest.df.cut),]

HR.KW.pest.df <- HR.pest.df[,c(as.character(HR.KW.nbrs), "OTU", "order", "family", "genus", "species")]
HR.KD.pest.df <- HR.pest.df[,c(as.character(HR.KD.nbrs), "OTU", "order", "family", "genus", "species")]
HR.AW.pest.df <- HR.pest.df[,c(as.character(HR.AW.nbrs), "OTU", "order", "family", "genus", "species")]
HR.AD.pest.df <- HR.pest.df[,c(as.character(HR.AD.nbrs), "OTU", "order", "family", "genus", "species")]
HR.BW.pest.df <- HR.pest.df[,c(as.character(HR.BW.nbrs), "OTU", "order", "family", "genus", "species")]
HR.BD.pest.df <- HR.pest.df[,c(as.character(HR.BD.nbrs), "OTU", "order", "family", "genus", "species")]

RA.KW.pest.df <- RA.pest.df[,c(as.character(RA.KW.nbrs), "OTU", "order", "family", "genus", "species")]
RA.KD.pest.df <- RA.pest.df[,c(as.character(RA.KD.nbrs), "OTU", "order", "family", "genus", "species")]
RA.AW.pest.df <- RA.pest.df[,c(as.character(RA.AW.nbrs), "OTU", "order", "family", "genus", "species")]
RA.AD.pest.df <- RA.pest.df[,c(as.character(RA.AD.nbrs), "OTU", "order", "family", "genus", "species")]
RA.BW.pest.df <- RA.pest.df[,c(as.character(RA.BW.nbrs), "OTU", "order", "family", "genus", "species")]
RA.BD.pest.df <- RA.pest.df[,c(as.character(RA.BD.nbrs), "OTU", "order", "family", "genus", "species")]


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
