#VCAP ----
load("/Volumes/Valeria/CRUK/AR_downstream_analysis/.RData")

# Probes correspondence
library(DBI)
library("illuminaHumanv3.db")
ExpAnnotation_VCap<-illuminaHumanv3fullReannotation()

annotation_vcap <- ExpAnnotation_VCap[which(ExpAnnotation_VCap$IlluminaID%in%probes$VCap.Probe_id),]
AnnotationVCap <- annotation_vcap[order(match(annotation_vcap$IlluminaID,probes$VCap.Probe_id)),]
rm(annotation_vcap)
rm(ExpAnnotation_VCap)

# AUTOCORRELATION ----------------------------------------------------------------------------------------

# Select the genes that have a behavioural change that differs from natural fluctuation.
# Autocorrelation higher than 0 (a threshold) indicate a change in behaviour possibly caused
# by AR 

#I calculate autocorrelation on each dataset and each repetition and then select the genes that pass the  
# autocorrelation filter in both repetitions

calculate_acf <- function(Data,ci){
  acf_data <- NULL
  for(i in 1:nrow(Data)){
    if(all(is.na(Data[i,]))){
      acf_data[i] <- NA
    }else{
      acf_data[i] <- acf(as.numeric(Data[i,]),type="correlation",lag.max=nrow(Data), plot=FALSE)$acf[2]
    }
  }
  # Only the positive acf gets selected. when acf<0 is because the signal goes up and down.
  acf_index <- which(acf_data>=ci) #which(abs(acf_data)>=ci) #
  acf_value <- acf_data[acf_index]
  acf_Data <- cbind(acf_index,acf_value)
  return(acf_Data)
}
# the qnorm... is how acf calculates the confidence intervals you can see in the plots.
# I use it as threshold.
VA1 <- calculate_acf(VCapAndrogen1,qnorm((1+0.95)/2)/sqrt(ncol(VCapAndrogen1)))
VA2 <- calculate_acf(VCapAndrogen2,qnorm((1+0.95)/2)/sqrt(ncol(VCapAndrogen2)))

VC1 <- calculate_acf(VCapControl1,qnorm((1+0.95)/2)/sqrt(ncol(VCapControl1)))
VC2 <- calculate_acf(VCapControl2,qnorm((1+0.95)/2)/sqrt(ncol(VCapControl2)))

# Select the genes that have been selected on both repetitions
index_VCapAndrogen <- intersect(VA1[,1],VA2[,1]) 

# Merge the indeces found in at least one repetition of Control.
index_VCapControl <- sort(union(VC1[,1],VC2[,1]))

# Compare these genes behaviour with the controls.
# Select them only if the gene is in the case but not in the Control.
index_VCap_Androgen_notControl <- sort(setdiff(index_VCapAndrogen,index_VCapControl))

# probes_acf_VCap <- probes[index_VCap_Androgen_notControl,]

# VARIANCE AND F-TEST - COMPARE VARIANCES OF TWO POPULATIONS --------------------------------------------

# Calculate variance for both androgens and controls. 
apply_ftest <- function(a1, a2, c1, c2, ind){
  var11 <- var12 <- var21 <- var22 <- NULL
  for(i in 1:nrow(a1)){
    var11[i] <- var.test(as.numeric(a1[i,]),as.numeric(c1[i,]))$p.value
    var12[i] <- var.test(as.numeric(a1[i,]),as.numeric(c2[i,]))$p.value
    var21[i] <- var.test(as.numeric(a2[i,]),as.numeric(c1[i,]))$p.value
    var22[i] <- var.test(as.numeric(a2[i,]),as.numeric(c2[i,]))$p.value  
  }
  
  #ind allow to refer the indeces to the original dataset
  a <- ind[which(var11<=0.05)]
  b <- ind[which(var12<=0.05)]
  d <- ind[which(var21<=0.05)]
  e <- ind[which(var22<=0.05)]
  
  ind_passed <- sort(unique(c(a,b,d,e)))
  return(ind_passed)
}

#indeces are still based on the original datasets
ind_acf_ftest_VCap <- apply_ftest(VCapAndrogen1[index_VCap_Androgen_notControl,], VCapAndrogen2[index_VCap_Androgen_notControl,],
                                  VCapControl1[index_VCap_Androgen_notControl,], VCapControl2[index_VCap_Androgen_notControl,],
                                  index_VCap_Androgen_notControl)

#Not for LNCap as this doesn't variates enough to measure the difference.

VCapAndrogen1_genesTF_passed <- VCapAndrogen1[ind_acf_ftest_VCap,]
VCapAndrogen2_genesTF_passed <- VCapAndrogen2[ind_acf_ftest_VCap,]
VCapControl1_genesTF_passed <- VCapControl1[ind_acf_ftest_VCap,]
VCapControl2_genesTF_passed <- VCapControl2[ind_acf_ftest_VCap,]

#from the probe file - contains ProbeType, VCapProbID, VCapSymbol, etc..
VCapProbes_genesTF_passed <- probes[ind_acf_ftest_VCap,c(1,2,3)]

#"hgnc_symbol-illumina_probe"
rownames(VCapAndrogen1_genesTF_passed) <- rownames(VCapAndrogen2_genesTF_passed) <- paste(VCapProbes_genesTF_passed[,2], VCapProbes_genesTF_passed[,3], sep=".")

#----------------------------------------------------------------------------------------------------------------------------------
# KEEP ALL GENES!!!! -------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

#Check if any of the genes selcted until now are transcription factors
# Script of this is in /media/bo01/     or    /Volumes/Valeria/CRUK/Human_TFs
# tfs_vcap19 <- get(load("/Volumes/Valeria/CRUK/Human_TFs/AnnotationVCap_TFs.RData"))[,c(1,11,12,13,16)]
annotation_vcap_genes_passed <- AnnotationVCap[which(AnnotationVCap$IlluminaID%in%VCapProbes_genesTF_passed$VCap.Probe_id),
                                               c(1,11,12,13,16)]
annotation_vcap_genes_passed_ordered <- annotation_vcap_genes_passed[order(match(annotation_vcap_genes_passed$IlluminaID,
                                                                                 VCapProbes_genesTF_passed$VCap.Probe_id)),]
library(tidyr)
annotation_vcap_multiple <- separate_rows(annotation_vcap_genes_passed_ordered, GenomicLocation, sep = ",")

# Convert into hg38 ----
#write the file
split_location <- function(Data){
  Split <- strsplit(as.character(Data$GenomicLocation),":")
  Chr <- sapply(Split, function(i) i[[1]])
  Chr <- gsub("chr", "", Chr)
  Start <- sapply(Split, function(i) i[[2]])
  End <- sapply(Split, function(i) i[[3]])
  Str <- sapply(Split, function(i) i[[4]])
  Split_location <- data.frame(EnsemblReannotated=Data$EnsemblReannotated,
                               chr=Chr,
                               start=as.integer(as.character(Start)),
                               end=as.integer(as.character(End)),
                               str=as.character(Str))
  return(Split_location)
}

tfs_vcap_lo <- split_location(annotation_vcap_multiple)
tfs_vcap_input_lo <- paste("chr",tfs_vcap_lo$chr,":",tfs_vcap_lo$start,"-",tfs_vcap_lo$end, sep="")
write.table(tfs_vcap_input_lo, file="genes_vcap_input_lo.bed",col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)

#Convert using liftover (UCSC)
# Go to the UCSC liftover load the input file and submit the query and then view 
# conversion will automatically download the file

#Read file
tfs_vcap_output_lo <- read.table("genes_vcap_output_lo.bed", header=FALSE)
tfs_vcap <- annotation_vcap_multiple
tfs_vcap$GenomicLocation <- as.character(tfs_vcap_output_lo[,1])
tfs_vcap$GenomicLocation <- gsub("-",":",tfs_vcap$GenomicLocation)
tfs_vcap$GenomicLocation <- paste(tfs_vcap$GenomicLocation,":",tfs_vcap_lo$str,sep="")

ind_tfs_VCap <- which(as.character(VCapProbes_genesTF_passed$VCap.Probe_id)%in%as.character(tfs_vcap$IlluminaID))

VCapAndrogen1_passed <- VCapAndrogen1_genesTF_passed[ind_tfs_VCap,]
VCapAndrogen2_passed <- VCapAndrogen2_genesTF_passed[ind_tfs_VCap,]
VCapControl1_passed <- VCapControl1_genesTF_passed[ind_tfs_VCap,]
VCapControl2_passed <- VCapControl2_genesTF_passed[ind_tfs_VCap,]

#from the probe file - contains ProbeType, VCapProbID, VCapSymbol, etc..
VCapProbes_passed <- VCapProbes_genesTF_passed[ind_tfs_VCap,]

#--------------------------------------------------------------------------------------------------------
# Ensembl_id_correspondence
tfs_vcap_in <- tfs_vcap[tfs_vcap$IlluminaID%in%as.character(VCapProbes_passed$VCap.Probe_id),]
ordered_tfs <- tfs_vcap_in[order(match(tfs_vcap_in$IlluminaID,VCapProbes_passed$VCap.Probe_id)),]
rm(tfs_vcap_in)

#--------------------------------------------------------------------------------------------------------
# TIME POINT OF CHANGE ----------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# Calculate where the behaviour of the gene change

#The repetitions have different time points therefore I cannot calculate the mean.
# I need to operate on each matrix separately.
# Calculate slope (y2-y1)/(x2-x1)

subtract_columns_matrix <- function(Matrix){  
  Matrix_sub <- matrix(0, nrow(Matrix), (ncol(Matrix)-1))  
  for(i in 1:(ncol(Matrix)-1)){  
    Matrix_sub[,i] <- Matrix[,i+1]-Matrix[,i]
  }
  return(Matrix_sub)  
}

subtract_elements_array <- function(Array){  
  Array_sub <- NULL  
  for(i in 1:(length(Array)-1)){  
    Array_sub[i] <- Array[i+1]-Array[i]
  }
  return(Array_sub)
}

divide_row_array <- function(Matrix, Array){
  Matrix_div <- matrix(0, nrow(Matrix), ncol(Matrix))    
  for(i in 1:nrow(Matrix)){
    Matrix_div[i,] <- round(Matrix[i,]/Array,digits=1)
  }
  return(Matrix_div)
}

mean_each2 <- function(Array){
  Array_mean <- NULL
  for(i in 1:(length(Array)-1)){
    Array_mean[i] <- round(mean(c(Array[i], Array[i+1])),digits=2)
  }
  return(Array_mean)
}

time_Androgen <- time_VCapAndrogen1
time_Control <- intersect(time_VCapControl1,time_VCapControl2)

# MEAN of the repetitions ----

VCapAndrogen_passed_mean <- VCapAndrogen1_passed 
VCapControl_passed_mean <- (VCapControl1_passed[,which(time_VCapControl1%in%time_Control)]+
                              VCapControl2_passed[,which(time_VCapControl2%in%time_Control)])/2

#LOESS of the mean of the repetitions -------------------------------------------------------------------
# I extrapolate the points of the controls using the time points of the androgen to not loose 
# all the time information

VCapAndrogen_passed_mean_loess <- VCapControl_passed_mean_loess <- matrix(0, 
                                                                          nrow(VCapAndrogen_passed_mean), 
                                                                          ncol(VCapAndrogen_passed_mean))
for(i in 1:nrow(VCapAndrogen_passed_mean)){
  VCapAndrogen_passed_mean_loess[i,] <- predict(loess(as.numeric(VCapAndrogen_passed_mean[i,])~time_Androgen, 
                                                      span=0.75))
  VCapControl_passed_mean_loess[i,] <- predict(loess(as.numeric(VCapControl_passed_mean[i,])~time_Control, 
                                                     span=0.75), time_Androgen)
}

# VCapAndrogen_passed_mean_loess_shifted <- cbind(VCapAndrogen_passed_mean_loess[,-1],
#                                                  VCapAndrogen_passed_mean_loess[,ncol(VCapAndrogen_passed_mean_loess)])
VCapAndrogen_passed_mean_shifted <- cbind(VCapAndrogen_passed_mean[,-1],
                                          VCapAndrogen_passed_mean[,ncol(VCapAndrogen_passed_mean)])

# Calculate the SLOPE
calculate_slope <- function(Matrix, Time){
  mat_sub <- subtract_columns_matrix(Matrix)
  time_sub <- subtract_elements_array(Time)
  matrix_slope <- divide_row_array(mat_sub, time_sub)
  return(matrix_slope)
}
VCapAndrogen_passed_mean_slope <- calculate_slope(VCapAndrogen_passed_mean_loess, time_Androgen)


#-------------------------------------------------------------------------------------------------------------------------------
library(TTR)

androgen_point_change <- point_change <- NULL
androgen_range_change <- range_change <- matrix(0,nrow(VCapAndrogen_passed_mean_loess),2)

thr_ccf <- round(qnorm((1+0.95)/2)/sqrt(ncol(VCapAndrogen_passed_mean_loess)), digits=3)

cc <- ac <- vector("list", nrow(VCapAndrogen_passed_mean_loess))
for(i in 1:nrow(VCapAndrogen_passed_mean_loess)){
  
  cc[[i]] <- ccf(as.numeric(VCapAndrogen_passed_mean_loess[i,]),
                 as.numeric(VCapControl_passed_mean_loess[i,]),
                 lag.max=length(as.numeric(VCapAndrogen_passed_mean_loess[i,])), 
                 type="correlation", plot=FALSE)
  
  # ac[[i]] <- acf(as.numeric(VCapAndrogen_passed_mean_loess[i,]),
  #           lag.max=length(as.numeric(VCapAndrogen_passed_mean_loess[i,])),  
  #           type="correlation", plot=FALSE)
  
  crosscorrelation <- round(cc[[i]]$acf,digits=3)
  
  # autocorrelation <- round(ac[[i]]$acf, digits=3)
  
  Lag <- cc[[i]]$lag
  
  above_thr <- which(abs(crosscorrelation) > 0)
  point_change_ccf <- abs(Lag[above_thr])
  
  thr_comp <- round(abs(VCapAndrogen_passed_mean_loess[i,1]-VCapControl_passed_mean_loess[i,1]),digits=1) + 0.2
  thr_ma <- which(round(abs(VCapAndrogen_passed_mean_loess[i,]-VCapControl_passed_mean_loess[i,]),digits=1)>thr_comp)
  mm <- round(abs(VCapAndrogen_passed_mean_loess[i,]-VCapControl_passed_mean_loess[i,]),digits=1)
  if(length(thr_ma)==0){
    thr_ma <- which(mm >= max(mm)-0.1)
  }
  point_change_ma <- thr_ma-1
  

  #NA if no time points intersect
  point_change[i] <- sort(intersect(point_change_ccf, point_change_ma))[1]
  
  range_change[i,1] <- ifelse(!is.na(point_change[i])&(point_change[i]==0), 
                                       point_change[i],
                                       point_change[i]-1)
  range_change[i,2] <- ifelse(!is.na(point_change[i])&(point_change[i]==0), 
                                       point_change[i]+2,
                                       point_change[i]+1)
  
  ##########

  maA <- round(SMA(as.numeric(VCapAndrogen_passed_mean_loess[i,]),3),1)
  ind_slope <- which(abs(VCapAndrogen_passed_mean_slope[i,])>0)
  ind_diff <- which(round(abs(maA - as.numeric(VCapAndrogen_passed_mean_loess[i,])),digits=1) >= 0.1)

  androgen_point_change[i] <- time_Androgen[intersect(ind_slope,ind_diff-2)[1]]
  
  androgen_range_change[i,1] <- ifelse(!is.na(androgen_point_change[i])&(androgen_point_change[i]==0), 
                                       androgen_point_change[i],
                                       androgen_point_change[i]-1)
  androgen_range_change[i,2] <- ifelse(!is.na(androgen_point_change[i])&(androgen_point_change[i]==0), 
                                       androgen_point_change[i]+2,
                                       androgen_point_change[i]+1)

}

dir.create("Figures")
pdf("Figures/compare_points_ranges_of_change.pdf",width = 7, height = 9)
par(mfrow=c(5,3),oma = c(2,2,2,2) + 0.1,mar = c(2,2,2,2) + 0.1)
for(i in 1:length(point_change)){
  Min <- min(c(as.numeric(VCapAndrogen_passed_mean[i,]),as.numeric(VCapControl_passed_mean[i,])))
  Max <- max(c(as.numeric(VCapAndrogen_passed_mean[i,]),as.numeric(VCapControl_passed_mean[i,])))
  plot(time_Androgen, VCapAndrogen_passed_mean[i,], 
       type="l", main=paste(i, VCapProbes_passed$VCap.Gene_symbol[i],sep="-"), 
       xlab=NULL, ylab=NULL, ylim=c(Min,Max))
  lines(time_Androgen, VCapAndrogen_passed_mean_loess[i,], col="blue")
  lines(time_Androgen, VCapControl_passed_mean_loess[i,], col="red")
  abline(v=point_change[i], col="green",lty=3)
  abline(v=range_change[i,1], col="green")
  abline(v=range_change[i,2], col="green")

  abline(v=androgen_point_change[i], col="pink",lty=3)
  abline(v=androgen_range_change[i,1], col="pink")
  abline(v=androgen_range_change[i,2], col="pink")
}
dev.off()


# Point change calculated using both androgen and control ----
time_change <- sort(unique(point_change))
# Indeces of the TFs that change at each time point
ind_change <- sapply(time_change, function(i) which(point_change==i), simplify=FALSE)
# Symbol_ProbeID of the TFs that change at each time point
tfs_time <- lapply(ind_change, function(i) rownames(VCapAndrogen_passed_mean)[i])

#Write the list of genes that change in a text file
tt <- 1
write.table(tfs_time[[tt]], file=paste("genes_that_change_at_", time_change[tt], ".txt", sep=""),
              quote=FALSE, col.names=FALSE, row.names=FALSE)

# Write a file with the list of genes and the corresponding point of change and androgen point of change ----

data_file <- cbind(VCapProbes_passed, point_change, androgen_point_change)
write.table(data_file, "genes_point_change_and_androgen_point.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Compare the Point change of the two approaches ----
difference_point_change <- point_change - androgen_point_change

hist(difference_point_change)
table(difference_point_change)


#----------------------------------------------------------------------------------------------------------------------------------
# Integrate CHIPSEQ DATA - apply_rcade --------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
# apply_rcade <- function(Data, Probes_passed, Probes_passed_names){
library(Rcade)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1- Calculate differential expression, fold change (microarray data) ----

library(limma)

Time <- c(time_VCapAndrogen1, time_VCapAndrogen2, time_VCapControl1, time_VCapControl2)
Group <- factor(c(rep("Androgen", length(c(time_VCapAndrogen1, time_VCapAndrogen2))),
                  rep("Control", length(c(time_VCapControl1, time_VCapControl2)))))
FileName <- c(paste("VA1_",time_VCapAndrogen1,sep=""),
              paste("VA2_",time_VCapAndrogen2,sep=""),
              paste("VC1_",time_VCapControl1,sep=""),
              paste("VC2_",time_VCapControl2,sep=""))

row.names(VCapAndrogen1_passed) <-  row.names(VCapAndrogen2_passed) <- as.character(VCapProbes_passed$VCap.Probe_id)
row.names(VCapControl1_passed) <-  row.names(VCapControl2_passed) <- as.character(VCapProbes_passed$VCap.Probe_id)
                                                                                    
                                                                                    
merged_VCap <- cbind(VCapAndrogen1_passed, VCapAndrogen2_passed, VCapControl1_passed, VCapControl2_passed)
colnames(merged_VCap) <- FileName

# both TFs and Genes ---

# merged_VCap <- cbind(VCapAndrogen1_genesTF_passed, VCapAndrogen2_genesTF_passed, VCapControl1_genesTF_passed, VCapControl2_genesTF_passed)
# colnames(merged_VCap) <- FileName

library(splines)
X <- ns(Time, df=5)

design <- model.matrix(~Group*X)

fit <- lmFit(merged_VCap,design)

fit <- eBayes(fit)

TopTable <- topTable(fit, adjust="BH", n=nrow(merged_VCap))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 3- Create the DE Matrix and the DE lookup table ----
tfs_unique <- unique(ordered_tfs[,c(1,5)])
gene_names <- tfs_unique[order(match(tfs_unique$IlluminaID,row.names(TopTable))),]

DE <- cbind(gene_names$EnsemblReannotated, TopTable)
colnames(DE) <- c("Ensembl.Gene.ID","logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
#Discard those genes with no corresponding Ensembl ID
DE <- DE[-which(is.na(DE$Ensembl.Gene.ID)),]

DElookup <- list(GeneID="Ensembl.Gene.ID", logFC="logFC", B="B")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4- Targets information - A matrix containing information about the .Bam les to be used in the analysis ----

samples_file <- read.csv("/Volumes/Valeria/CRUK/AR_all_folder/AR_final_code/Dora_AR/samples.csv", 
                         header=TRUE, colClasses = c(NA,"NULL","NULL","NULL","NULL","NULL",NA,NA,NA,"NULL","NULL",NA,NA)) 

fileID <- as.character(samples_file$SRA)
sampleID <- samples_file$Experiment
factorID <- samples_file$Sample.Input
levels(factorID) <- c("Input", "S")

# 4a- Order the files' name in the same order as in the table ----
bam_bai <- dir("/Volumes/Valeria/CRUK/AR_all_folder/AR_final_code/Dora_AR/bam", pattern = ".bam")
file_dir <- bam_bai[-c(grep("bai", bam_bai),grep("samstat", bam_bai))]

file_dir_name <- sapply(strsplit(file_dir,"_"), function(i) i[1])

#x nell'ordine di y --> x[order(match(x,y))]
filepath <- file_dir[order(match(file_dir_name,fileID))]

targets <- data.frame(fileID,sampleID,factorID,filepath)

# 4b- Use only the samples you are interested in analysing ----
interest <- c("SRR039769", "SRR039774", "SRR039773", "SRR039775")
targets_int <- targets[which(targets$fileID%in%interest),]
# write.table(targets_int, file="/Volumes/Valeria/CRUK/AR_downstream_analysis/EncodeData/targets_chipseq_peaks_id.txt", 
#             quote=FALSE, row.names=FALSE, col.names=FALSE)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5- Annotation information ----

ordered_tfs_position2 <- split_location(ordered_tfs)

ordered_tfs_position <- ordered_tfs_position2[-which(is.na(ordered_tfs_position2[,1])),]
colnames(ordered_tfs_position) <- c("Ensembl.Gene.ID","chr","start", "end", "str")

# a <- ordered_tfs_position[20:25,]
ChIPannoZones_1500 <- defineBins(ordered_tfs_position, zone=c(-1500, 1500), geneID="Ensembl.Gene.ID")
ChIPannoZones_5 <- defineBins(ordered_tfs_position, zone=c(-5000, 5000), geneID="Ensembl.Gene.ID")
ChIPannoZones_10 <- defineBins(ordered_tfs_position, zone=c(-10000, 10000), geneID="Ensembl.Gene.ID")
ChIPannoZones_25 <- defineBins(ordered_tfs_position, zone=c(-25000, 25000), geneID="Ensembl.Gene.ID")

# levels(ordered_tfs_position2$str) <- c("-1","1")
# ordered_tfs_position2$str <- as.character(ordered_tfs_position2$str)
# ordered_tfs_position2$str[which(ordered_tfs_position2$str=="+")] <- 1
# ordered_tfs_position2$str[which(ordered_tfs_position2$str=="-")] <- -1
# 
# ordered_tfs_position <- ordered_tfs_position2[-which(is.na(ordered_tfs_position2[,1])),]
# colnames(ordered_tfs_position) <- c("Ensembl.Gene.ID","chr","start", "end", "str")
# 
# ind_file <- 1:nrow(ordered_tfs_position)
# file_to_write <- data.frame(ordered_tfs_position[,c(2,3,4)],ind_file)
# write.table(file_to_write, 
#             file="/Volumes/Valeria/CRUK/AR_downstream_analysis/EncodeData/genes_location_input_RCade.bed", 
#             sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 6- Prior specification ----
DE.prior <- 0.01
prior.mode <- "keepChIP"
prior <- c("D|C"=0.05, "D|notC"=0.005)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 7- Analysis function ----
Dir <- "/Volumes/Valeria/CRUK/AR_all_folder/AR_final_code/Dora_AR/bam"
cl <- NULL
dir.create("RCade")

run_rcade <- function(ChIPannoZones, Number){
  Rcade <- RcadeAnalysis(DE, ChIPannoZones, annoZoneGeneidName = "Ensembl.Gene.ID", 
                         ChIPtargets = targets_int, ChIPfileDir = Dir, cl=cl, 
                         DE.prior=DE.prior, prior.mode=prior.mode, prior=prior, DElookup=DElookup)
  dir.create(paste("RCade/Rcade_files_B_",Number,sep=""))
  exportRcade(Rcade, directory = paste("RCade/Rcade_files_B_",Number,sep=""),  cutoffMode = "B", cutoffArg = -5, removeDuplicates = "beforeCutoff")
  
  # dir.create(paste("RCade/Rcade_files_FDR_",Number,sep=""))
  # exportRcade(Rcade, directory = paste("RCade/Rcade_files_FDR_",Number,sep=""),  cutoffMode = "FDR", removeDuplicates = "beforeCutoff")
  
  rcade_results <- read.csv(paste("RCade/Rcade_files_B_",Number,"/DEandChIP.csv",sep=""), header=TRUE, colClasses = c(NA,rep("NULL",11),NA)) 
  # rcade_output <- getRcade(Rcade)
  return(rcade_results)
}

#Genes within the 1500 bases window ----
rcade_results_1500 <- run_rcade(ChIPannoZones_1500,"1500")
dim(rcade_results_1500)

symbol_1500 <- annotation_vcap_multiple[which(annotation_vcap_multiple$EnsemblReannotated%in%rcade_results_1500$geneID),c("SymbolReannotated","EnsemblReannotated")]
symbol_1500_ordered <- unique(symbol_1500[order(match(symbol_1500$EnsemblReannotated,rcade_results_1500$geneID)),])
dim(symbol_1500_ordered)
write.table(symbol_1500_ordered, "symbol_1500_ordered.txt",quote = FALSE, row.names = FALSE)

#Genes within the 5000 bases window ----
rcade_results_5 <- run_rcade(ChIPannoZones_5,"5")
dim(rcade_results_5)

symbol_5000 <- annotation_vcap_multiple[which(annotation_vcap_multiple$EnsemblReannotated%in%rcade_results_5$geneID),c("SymbolReannotated","EnsemblReannotated")]
symbol_5000_ordered <- unique(symbol_5000[order(match(symbol_5000$EnsemblReannotated,rcade_results_5$geneID)),])
dim(symbol_5000_ordered)
write.table(symbol_5000_ordered, "symbol_5000_ordered.txt",quote = FALSE, row.names = FALSE)

#Genes within the 10000 bases window ----
rcade_results_10 <- run_rcade(ChIPannoZones_10,"10")
dim(rcade_results_10)

symbol_10000 <- annotation_vcap_multiple[which(annotation_vcap_multiple$EnsemblReannotated%in%rcade_results_10$geneID),c("SymbolReannotated","EnsemblReannotated")]
symbol_10000_ordered <- unique(symbol_10000[order(match(symbol_10000$EnsemblReannotated,rcade_results_10$geneID)),])
dim(symbol_10000_ordered)
write.table(symbol_10000_ordered, "symbol_10000_ordered.txt",quote = FALSE, row.names = FALSE)

#Genes within the 25000 bases window ----
rcade_results_25 <- run_rcade(ChIPannoZones_25,"25")
dim(rcade_results_25)

symbol_25000 <- annotation_vcap_multiple[which(annotation_vcap_multiple$EnsemblReannotated%in%rcade_results_25$geneID),c("SymbolReannotated","EnsemblReannotated")]
symbol_25000_ordered <- unique(symbol_25000[order(match(symbol_25000$EnsemblReannotated,rcade_results_25$geneID)),])
dim(symbol_25000_ordered)
write.table(symbol_25000_ordered, "symbol_25000_ordered.txt",quote = FALSE, row.names = FALSE)


genes_in_all_windows <- intersect(symbol_25000_ordered$SymbolReannotated,intersect(symbol_10000_ordered$SymbolReannotated,intersect(symbol_5000_ordered$SymbolReannotated,symbol_1500_ordered$SymbolReannotated)))
write.table(genes_in_all_windows, "symbol_genes_in_all_windows.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

#Compare 1.5 kb and 5kb ----
# Genes in 1500 but not in 5000
genes_in1500_not5000 <- symbol_1500_ordered$SymbolReannotated[which(!symbol_1500_ordered$SymbolReannotated%in%symbol_5000_ordered$SymbolReannotated)]
length(genes_in1500_not5000)
# Genes in 5000 but not in 1500
genes_in5000_not1500 <- symbol_5000_ordered$SymbolReannotated[which(!symbol_5000_ordered$SymbolReannotated%in%symbol_1500_ordered$SymbolReannotated)]
length(genes_in5000_not1500)
# Genes in both 1.5kb and 5kb
genes_in5000_and1500 <- symbol_5000_ordered$SymbolReannotated[which(symbol_5000_ordered$SymbolReannotated%in%symbol_1500_ordered$SymbolReannotated)]
length(genes_in5000_and1500)

#Compare 1.5 kb and 10kb ----
# Genes in 5000 but not in 10000
genes_in1500_not10000 <- symbol_1500_ordered$SymbolReannotated[which(!symbol_1500_ordered$SymbolReannotated%in%symbol_10000_ordered$SymbolReannotated)]
length(genes_in1500_not10000)
# Genes in 10000 but not in 1500
genes_in10000_not1500 <- symbol_10000_ordered$SymbolReannotated[which(!symbol_10000_ordered$SymbolReannotated%in%symbol_1500_ordered$SymbolReannotated)]
length(genes_in10000_not1500)
# Genes in both 1.5kb and 10kb
genes_in10000_and1500 <- symbol_10000_ordered$SymbolReannotated[which(symbol_10000_ordered$SymbolReannotated%in%symbol_1500_ordered$SymbolReannotated)]
length(genes_in10000_and1500)

#Compare 1.5 kb and 25kb ----
# Genes in 1500 but not in 25000
genes_in1500_not25000 <- symbol_1500_ordered$SymbolReannotated[which(!symbol_1500_ordered$SymbolReannotated%in%symbol_25000_ordered$SymbolReannotated)]
write.table(genes_in1500_not25000, "genes_in1500_not25000.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
# Genes in 25000 but not in 1500
genes_in25000_not1500 <- symbol_25000_ordered$SymbolReannotated[which(!symbol_25000_ordered$SymbolReannotated%in%symbol_1500_ordered$SymbolReannotated)]
write.table(genes_in25000_not1500, "genes_in25000_not1500.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
# Genes in both 1.5kb and 25kb
genes_in25000_and1500 <- symbol_25000_ordered$SymbolReannotated[which(symbol_25000_ordered$SymbolReannotated%in%symbol_1500_ordered$SymbolReannotated)]
write.table(genes_in25000_and1500, "genes_in25000_and1500.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

length(genes_in25000_and1500)/length(unique(symbol_25000_ordered$SymbolReannotated))
length(genes_in25000_and1500)/length(unique(symbol_1500_ordered$SymbolReannotated))

