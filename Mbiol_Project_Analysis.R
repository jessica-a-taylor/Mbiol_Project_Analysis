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
# seq.data <- readxl::read_xlsx(choose.files())

# filter for insect OTUs
# insect.data <- seq.data[seq.data$class=="Insecta",]
# insect.data <- insect.data[,-c(56,185,272,274:291)]

# removing OTUs with less than 10 reads and those represented by less than 
#     1% of the reads for each sample

#####
#   insect.data.cut <- c()                                          
#   for (row in 1:nrow(insect.data)) {                              
#     D <- FALSE                                                    
#     for (col in 1:ncol(insect.data[,2:270])) {                    
#       ifelse(insect.data[row, col] >= 10,                         
#              D <- TRUE, D <- D)                                   
#     } 
#     if (D == FALSE) {                                             
#       insect.data.cut <- append(insect.data.cut, row)             
#     }                                                             
#   }                                                               
#                                                                   
#   insect.data <- insect.data[-c(insect.data.cut),]                
#                                                                   
#   new.insect.data <- insect.data[,-c(1,271:276)] 
#
#   sum.reads <- c()                                                
#   for (nid.sum in 1:ncol(new.insect.data)) {                      
#     sum.reads <- append(sum.reads, sum(new.insect.data[,nid.sum]))
#   }                                                               
#                                                                   
#   new.insect.data <- new.insect.data[,-c(247,248)]
#   sum.reads <- c()                                                
#   for (nid.sum in 1:ncol(new.insect.data)) {                      
#     sum.reads <- append(sum.reads, sum(new.insect.data[,nid.sum]))
#   }                                                               
#                                                                   
#   colnames(new.insect.data) <- c(sum.reads) 
#                                                                   
#   for (idc in as.character(sum.reads)) {                          
#     normalised.diet <- (new.insect.data[,idc]/as.numeric(idc))*100
#     new.insect.data[,idc] <- normalised.diet                      
#   }                                                               
#                                                                   
#   new.insect.data.cut <- c()                                      
#   for (row in 1:nrow(new.insect.data)) {
#     D <- FALSE                                                    
#     for (col in 1:ncol(new.insect.data)) {                        
#       ifelse(new.insect.data[row, col] >= 1,                      
#              D <- TRUE, D <- D)                                   
#     }                                                             
#     if (D == FALSE) {                                             
#       new.insect.data.cut <- append(new.insect.data.cut, row)     
#     }                                                             
#   }

#   insect.data <- insect.data[-c(new.insect.data.cut),]            
#   insect.data <- insect.data[,-c(248,249)] 
#####

# diet calculations
# new.insect.data <- as.data.frame(insect.data[,-c(1,269:274)])

# function for calculating FOO
#####
#   FOO.function <- function() {                                                                
#                                                                                               
#     insect.FOO.counts <- c()                                                                  
#     for (seq_row in 1:nrow(new.insect.data)) {                                                
#       I <- 0                                                                                  
#       for (seq_col in 1:ncol(new.insect.data)) {                                              
#         ifelse(new.insect.data[seq_row, seq_col] > 0,                                         
#                I <- I + 1, I <- I) 
#       }                                                                                       
#       insect.FOO.counts <- append(insect.FOO.counts, (1/length(names(new.insect.data)))*I*100)
#     }                                                                                         
#   }
#####

# FOO.function()
# insect.data <- cbind(insect.data, insect.FOO.counts)

# function for calculating POO
#####
#   POO.function <- function() {                                      
#     insect.POO.counts <- c()                                        
#     for (seq_row in 1:nrow(new.insect.data)) {                      
#       I <- 0                                                        
#       for (seq_col in 1:ncol(new.insect.data)) {                    
#         ifelse(new.insect.data[seq_row, seq_col] > 0,               
#                I <- I + 1, I <- I)                                  
#       }  
#                                                                     
#     sum.insect.POO.counts <- sum(insect.POO.counts)                 
#     insect.POO <- c()                                               
#     for (i in insect.POO.counts) {                                  
#       insect.POO <- append(insect.POO,(i/sum.insect.POO.counts)*100)
#     }                                                               
#   }
#####

# POO.function()
# insect.data <- cbind(insect.data, insect.POO)

# function for calculating wPOO
#####
#   wPOO.function <- function() {                                                                              
#     for (row.no in 1:nrow(new.insect.data)) {                                                                
#       for (col.no in 1:ncol(new.insect.data)) {                                                              
#         if (new.insect.data[row.no, col.no] >0) {                                                            
#           new.insect.data[row.no, col.no] <- 1
#         }                                                                                                    
#       }                                                                                                      
#     
#                                                                                                              
#     insect.sample.diet <- c()                                                                                
#     for (seq_col in 1:ncol(new.insect.data)) {                                                               
#       J <- 0
#       for (seq_row in 1:nrow(new.insect.data)) {                                                             
#         ifelse(new.insect.data[seq_row, seq_col] > 0, J <- J + 1, J <- J)                                    
#       }                                                                                                      
#       insect.sample.diet <- append(insect.sample.diet, J)                                                    
#     }
#     colnames(new.insect.data) <- c(insect.sample.diet)                                                       
#     insect.value.df <- new.insect.data                                                                       
#                                                                                                              
#     new.insect.data <- as.data.frame(new.insect.data)                                                        
#     reads <- c()
#     for (insect.col in 1:ncol(new.insect.data)) {                                                            
#       reads <- append(reads, new.insect.data[,insect.col])                                                   
#     }                                                                                                        
#
#     insect.df <- data.frame(Row = rep(1:950, times = 267),                                                   
#                             Reads = reads,                                                                   
#                             Total.reads = rep(insect.sample.diet, each = 950))                               
#
#     insect.df <- cbind(insect.df, insect.df$Reads/insect.df$Total.reads)                                     
#                                                                                                              
#     insect.value.df <- new.insect.data                                                                       
#     for (row in 1:nrow(insect.value.df)) {
#       new.insect.df <- insect.df[insect.df$Row==row,]                                                        
#       insect.value.df[row,] <- new.insect.df$`insect.df$Reads/insect.df$Total.reads`                         
#     }                                                                                                        
#
#     insect.wPOO <- c()                                                                                       
#     for (value_row in 1:nrow(insect.value.df)) {                                                             
#       insect.wPOO <- append(insect.wPOO, 100*(1/length(insect.sample.diet))*sum(insect.value.df[value_row,]))
#     }                                                                                                        
#   }
#####

# wPOO.function()
# insect.data <- cbind(insect.data, insect.wPOO)

# writexl::write_xlsx(insect.data, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\insect.filtered.data.xlsx")

# import habitat data
# select columns of interest and remove NAs
#####
#   zbj_data <- readxl::read_xlsx(choose.files())                                        
#   farms <- readxl::read_xlsx(choose.files(), sheet = "Samples")                        
#   farms <- data.frame(farms$Lab.nbr...26, farms$animal, farms$Species, farms$Location, 
#                       farms$Site, as.character.POSIXt(farms$Date))                     
#   farms <- farms[-145,]                                                                
#                                                                                        
#   farms.na <- which(is.na(farms$farms.Lab.nbr...26), arr.ind=TRUE)                     
#   farms <- farms[-c(farms.na),]
#                                                                                        
#   sort_farms <- farms[c(4,5)]                                                          
#   farm_names <- unique(sort_farms$farms.Location)                                      
#   farm_names <- farm_names[!is.na(farm_names)]
#####

# use hash match up missing farm names
#####
#   nameMap <- hash()                                
#                                                    
#   for (r in 1:nrow(sort_farms)) {                  
#     locationName <- sort_farms[r, "farms.Location"]
#     siteName <- sort_farms[r, "farms.Site"]        
#     nameMap[[siteName]] <- locationName            
#   }
#####

# find lab numbers common to both datasets to identify which samples were actually sequenced
#####
#   farms.ln <- as.character(farms$farms.Lab.nbr...26)
#                                                     
#   zbj_names <- names(insect.data[,2:268])           
#                                                     
#   Lab.no <- c()                                     
#   for (n in zbj_names) {                            
#     for (l in farms.ln) {                           
#       if (l==n)
#         Lab.no <- append(Lab.no, as.numeric(n))     
#     }                                               
#   }
#####

# compile information of interest from both datasets
#####
#   zbj_df <- data.frame()                                                                                                                  
#                                                                                                                                           
#   for (i in Lab.no) {                                                                                                                     
#     seqSamples <- farms[which(farms$farms.Lab.nbr...26==i),]                                                                              
#     zbj_df <- rbind(zbj_df, seqSamples) 
#   }                                                                                                                                       
#   for (d in 1:nrow(zbj_df)) {                                                                                                             
#     if (is.na(zbj_df[d, "farms.Location"])) {                                                                                             
#       zbj_df[d, "farms.Location"] <- nameMap[[zbj_df[d, "farms.Site"]]]                                                                   
#     }
#     zbj_df[d, "as.character.POSIXt.farms.Date."] <- ifelse(str_detect(zbj_df[d, "as.character.POSIXt.farms.Date."], "2017"), "wet", "dry")
#   }
#                                                                                                                                           
#   zbj_final <- zbj_df[order(zbj_df[,4]),]                                                                                                 
#   colnames(zbj_final) <- c("Lab.nbr", "Animal", "Species", "Location", "Site", "Season")
#####

# write_xlsx(zbj_final, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\zbj_final.xlsx")

insect.data <- readxl::read_xlsx("C:\Users\jexy2\OneDrive\Documents\Mbiol project\insect.filtered.data.xlsx")
zbj_final <- readxl::read_xlsx("C:\Users\jexy2\OneDrive\Documents\Mbiol project\zbj_final.xlsx")

# create hash for matching lab numbers to species name
#####
#   predator.names <- hash()                             
#   for (row in 1:nrow(zbj_final)) {                     
#     pred.nbr <- as.character(zbj_final[row, "Lab.nbr"])
#     pred.sp <- as.character(zbj_final[row,"Species"])  
#     predator.names[[pred.nbr]] <- pred.sp              
#   }
#####

# use hash to create a matrix of sequence data for the focal species
#####
#   zbj.insect.data <- insect.data[,c(2:268)]                                                  
#   zbj.insect.data <- zbj.insect.data[,-c(195, 218, 222)]                                     
#   new.zbj.insect.data <- zbj.insect.data                                                     
#                                                                                              
#   pred.colnames <- c()                                                                       
#   for (name in names(new.zbj.insect.data)) {                                                 
#     pred.colnames <- append(pred.colnames, predator.names[[name]])                           
#   }
#
#   colnames(new.zbj.insect.data) <- c(pred.colnames)                                          
#                                                                                              
#   HR.focal.nbrs <- which(names(new.zbj.insect.data) == "Hipposideros ruber", arr.ind = TRUE) 
#   RA.focal.nbrs <- which(names(new.zbj.insect.data) == "Rhinolophus alcyone", arr.ind = TRUE)
#                                                                                              
#   all.focal.nbrs <- c(HR.focal.nbrs, RA.focal.nbrs)                                          
#   all.focal.nbrs <- sort(all.focal.nbrs)                                                     
#                                                                                              
#   zbj.insect.data.names <- names(zbj.insect.data)
#   focal.nbrs <- zbj.insect.data.names[all.focal.nbrs]                                        
#   zbj.focal <- zbj.insect.data[,c(focal.nbrs)]                                               
#####

# flip columns and rows
# zbj.focal <- t(zbj.focal)

# convert seq data to occurrence data
#####
#   zbj.focal <- t(zbj.focal)         
#                                     
#   for (row in 1:nrow(zbj.focal)) {  
#     for (col in 1:ncol(zbj.focal)) {
#       if(zbj.focal[row,col] > 0) {  
#         zbj.focal[row,col] <- 1     
#       }                             
#     }                               
#   }
#####

# add columns with sample and habitat info
#####                                                                                          
#   new.zbj.insect.data.names <- names(new.zbj.insect.data)                                    
#   focal.names <- new.zbj.insect.data.names[all.focal.nbrs]
#   zbj.focal <- as.data.frame(zbj.focal)                       
#   zbj.focal <- cbind(zbj.focal, focal.nbrs)                   
#   zbj.focal <- cbind(zbj.focal, focal.names)                  
#                                                               
#   focal.info <- data.frame()                                  
#   for (nbr in as.numeric(focal.nbrs)) {                       
#     for (row in 1:nrow(zbj_final)) {                          
#       ifelse(zbj_final[row, "Lab.nbr"] == nbr,                
#              focal.info <- rbind(focal.info, zbj_final[row,]),
#              focal.info <- focal.info)                        
#     }                                                         
#   }                                                           
#                                                               
#   zbj.focal <- cbind(zbj.focal, focal.info[,c(4,6)])
#   focal.info <- focal.info[,-c(1,2,5)]
#   colnames(zbj.focal) <- c(insect.data$OTU, "Lab.nbr", "Species", "Location", "Season")
#####

# write_xlsx(zbj.focal, "C:\\Users\\jexy2\\OneDrive\\Documents\\Mbiol project\\focal.species.data.xlsx")

# db-rda analysis
#####
#   library(vegan)                                                                              
#   library(dplyr)                                                                              
#   library(ggplot2)                                                                            
#   library(ggalt)                                                                              
#   library(ggforce)                                                                            
#   library(concaveman)                                                                         
#                                                                                               
#   focal.nmds <- metaMDS(zbj.focal[1:950], distance = "jaccard")
#   plot(focal.nmds, type = "t", main = paste("NMDS/Jaccard - Stress = ",                       
#                                             round(focal.nmds$stress, 3)))                     
#                                                                                               
#   focal.envfit <- envfit(focal.nmds, focal.info, permutations = 999, na.rm = TRUE)            
#   plot(focal.envfit)                                                                          
#                                                                                               
#   en_coord_cat <- as.data.frame(scores(focal.envfit, "factors")) #* ordiArrowMul(focal.envfit)
#                                                                                               
#   data.scores <- as.data.frame(scores(focal.nmds))
#   data.scores <- cbind(data.scores, focal.info)                                               
#                                                                                               
#   data.scores <- data.scores[-c(155, 158),]
#####

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

focal.interact.species <- capscale(zbj.focal[1:950] ~ zbj.focal$Species, distance = "jaccard", add = TRUE)
anova(focal.interact.species, step = 1000, perm.max = 1000)

focal.interact.season <- capscale(zbj.focal[1:950] ~ zbj.focal$Season, distance = "jaccard", add = TRUE)
anova(focal.interact.season, step = 1000, perm.max = 1000)

focal.interact.location <- capscale(zbj.focal[1:950] ~ zbj.focal$Location, distance = "jaccard", add = TRUE)
anova(focal.interact.location, step = 1000, perm.max = 1000)

focal.interact.all <- capscale(zbj.focal[1:950] ~ zbj.focal$Species + zbj.focal$Season + zbj.focal$Location, distance = "jaccard", add = TRUE)
plot(focal.interact.all)
anova(focal.interact.all, step = 1000, perm.max = 1000)
anova(focal.interact.all, by="axis", perm.max=500)
anova(focal.interact.all, by="terms", permu=200)
