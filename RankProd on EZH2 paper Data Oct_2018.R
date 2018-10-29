###############################################################################################
# Demonstration of RankProd via re-analysis of EZH2 paper data (10.18632/oncotarget.25128)
###############################################################################################


#experimental design and shorthand key:
  #Two cell lines: MMM1 ("M") and FLAM76 ("F")
  #HMCLs harvested on three time points (days: "1", "4" and "5")
  #HMCLs were treated with:
    #media/ctrl ("A"; days 1,4 & 5)
    #EPZ6438 (EZH2i) ("E"; days 1,4 & 5)
    #Panobinostat (pan-HDACi) ("B"; day 5)
    #combination ("G"; day 5)
  #Triplicate samples ("1","2" & "3")
  #example: M_A4_2 (cellline_treatment/day_replicate)


#RNA-seq data presented in the above manuscript quantified, merged and calculated differential expression in one step 
#To generate quantified FPKM values for each replicate individually: BAM files (aligned to hg19 via tophat) were re-quantified with Cufflinks v2.2.1 via UMN Galaxy
#cufflinks was chosen for consistancy towards comparing RankProd and cuffdiff outputs downstream
#output files are tab-delimited - manually downloaded from galaxy history


#generate a annotation Key for cufflinks output files (n=48)
  #"RQD" = re-quantified data (individual replicates as opposed to merged cuffdiff output)
  RQD_AnnotationNames <- character(length = 48)
  a<-1
  for(b in c("F","M")){
    for (c in c("A1","A4","A5","B5","E1","E4","E5","G5")) {
      for(d in c(1,2,3)){
        RQD_AnnotationNames[a]<-unlist(paste(b,c,d,sep = "_"))
        a<-a+1
      }
    }
  }
  rm(a,b,c,d)
  
  RQD_AnnotationKey <- data.frame(sample.name = RQD_AnnotationNames, 
                                  cell.line = rep(c("F","M"), each = 24),
                                  treatment = rep(c("A","A","A","B","E","E","E","G"),each = 3, times = 2),
                                  day = rep(c("1","4","5","5"), each = 3, times = 4),
                                  replicate = rep(c("1","2","3"),16)
                                  )
                                  

#parse in .tabular files output by cufflinks
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data cufflinks output")
  RQD_filenames <- list.files() #filenames are "01.tabular"..."48.tabular"
  for (x in c(1:48)) {
    assign(RQD_AnnotationNames[x],read.table(RQD_filenames[x], sep = '\t',header = TRUE, stringsAsFactors = FALSE))
  }
  rm(x)

#resolve duplicate gene_id's (known bug in cufflinks) and simplify data.frame
  #source: separate transcript clusters map to the same gene
  #solution: sum FPKMs (assumption: different clusters are not biologically distinct)
  RQD_unique.geneIDs <- sort(unique(get(RQD_AnnotationNames[1])$gene_id)) #gene names are shared in all samples but are not ordered identically
  for (x in RQD_AnnotationNames) {
    print(x)
    df.in <- get(x)
    df.out <- data.frame(row.names = RQD_unique.geneIDs)
    df.out$FPKM <- 0
    for (y in c(1:length(RQD_unique.geneIDs))) {
      df.out$FPKM[y]<- sum(df.in$FPKM[df.in$gene_id == RQD_unique.geneIDs[y]])
    }
    assign(x,df.out)
  }
  rm(x,y,df.in,df.out)
  
#filtering
  #remove genes with less than 1 FPKM in every sample
  filter <- rep(0,23288)
  for (x in RQD_AnnotationNames) {
     filter <- filter + (get(x)$FPKM >= 1) # querries all samples, where filter remains equal to 0 that gene is <1 FPKM in all samples
  }
  rm(x)
  table(filter) #view distribution - 8888 genes are <1 FPKM in every sample (8579 genes have FPKM >= 1 in every sample)
  for (x in RQD_AnnotationNames) {
    filtered_df <- get(x)
    filtered_df <- data.frame(row.names = RQD_unique.geneIDs[(filter != 0)], FPKM = filtered_df[(filter != 0),])
    assign(x,filtered_df)
  }
  RQD_filtered.geneIDs <- RQD_unique.geneIDs[(filter != 0)]
  rm(x,filtered_df)
  
############################################
#RankProd 2.0
  
#publication (10.1093/bioinformatics/btx292)
#documentation & tutorial (https://www.bioconductor.org/packages/devel/bioc/html/RankProd.html)
  
#install
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("RankProd")

library(RankProd)

#note: these data are consistant with 'single origin' data as described in documentation

#define function to generate argument for RankProducts function 
  #returns FPKM data.frame with samples as cols, genes as rows
  #expects shorthand outlined above (as char)
  #'treatment1' and 'treatment2' are two treatments being compared (assume same day)
  RP.expr.data <- function(cell.line.out, treatment1.out,treatment2.out, day.out){ 
    names1 <- as.character(subset(RQD_AnnotationKey,cell.line == cell.line.out & treatment == treatment1.out & day == day.out)$sample.name)
    names2 <- as.character(subset(RQD_AnnotationKey,cell.line == cell.line.out & treatment == treatment2.out & day == day.out)$sample.name)
    names.out <- c(names1,names2)
    print(names1)
    data.out <- data.frame(row.names = row.names(get(names.out[1])))
    for(x in names.out){
      data.out[,x]<-get(x)$FPKM
    }
  return(data.out)  
  }
  

#cl argument will be the same in all these cases
    #vector separating the two test group 0 = ctrl, 1 = 'treated'
    #length must equal ncol of data
cl <- c(0,0,0,1,1,1) 

#differential expression via RankProducts function
  #assumed log transformed expr values by default - these data ARE NOT!
  #imputes NA values by default (no NAs in these data)
  #calculates rank product (vs rank sum) by default
  #try on MA5 vs MG5  
  RP.out <- RankProducts(RP.expr.data("M","A","G","5"),cl,logged = FALSE,gene.names = RQD_filtered.geneIDs,rand=123) #rand value set for reproducability for now
  str(RP.out) #output structured as list of lists

#plot results
  plotRP(RP.out, cutoff=0.05)
  
#print sigdiff genes
  #note: pfp cutoff method much more conservitave than pval without p.adjust (in MA5vMG5 produces 322 and 3325 sigdiff genes respectively - becomes 328 with BH adjust)
  topGene(RP.out,cutoff = 0.05, method = "pfp",logged = FALSE, gene.names = RQD_filtered.geneIDs)

#systematically evaluate many diffexpr comparisons and creat custom output comparing to cuffdiff output used in paper
  #declare func to generate merged output as data.frame
  RP.custom <- function(cell.line.in, treatment1.in, treatment2.in, day.in,FPKM.filtercuttoff = 0,FC.filtercuttoff = 1.5){
    RP.custom.out <- data.frame(row.names = RQD_filtered.geneIDs)
    RP.custom.in <- RankProducts(RP.expr.data(cell.line.in, treatment1.in, treatment2.in, day.in),cl,logged = FALSE,gene.names = RQD_filtered.geneIDs) #random seed not fixed as above
    print(topGene(RP.custom.in,cutoff = 0.05, method = "pfp",logged = FALSE, gene.names = RQD_filtered.geneIDs))
    print("pfp<0.05")
    print(sum(RP.custom.in$pfp <0.05))
    #(cuffdiff outputs previously parsed into workspace labeled with following format: e.g. MA1vME1 or FA5vFB5 )
    #see below to parse in
    Old.data.name <- paste(cell.line.in,treatment1.in,day.in,"v",cell.line.in,treatment2.in,day.in, sep = "")
    print(Old.data.name)
    Cuffdiff.out <- get(Old.data.name) #these data.frames do not have duplicate gene names
    #all gene names perfectly between cufflinks/cuffdiff outputs > sum(RQD_filtered.geneIDs %in% MA5vMB5$test_id)
    #subset cuffdiff output to match querried genes
    row.names(Cuffdiff.out)<-Cuffdiff.out$gene
    Cuffdiff.out <- Cuffdiff.out[RQD_filtered.geneIDs,]
    #assume genes cant be both sigdiff up and down via RP
    RP.updownkey <- data.frame ("up" = (1/RP.custom.in$AveFC)>1,"down" = (1/RP.custom.in$AveFC)<1)
    RP.custom.out$RP.pfp <- rowSums(RP.updownkey*RP.custom.in$pfp)
    RP.custom.out$RP.pval <- rowSums(RP.updownkey*RP.custom.in$pval)
    RP.custom.out$RP.qval <- rowSums(RP.updownkey*p.adjust(RP.custom.in$pval, method = "BH"))
    RP.custom.out$RP.fc <- 1/RP.custom.in$AveFC # invert so that >1 = upreg and <1 = downreg (class2/class1) as in cuffdiff output
    RP.custom.out$CD.pval <- Cuffdiff.out$p_value
    RP.custom.out$CD.qval <- Cuffdiff.out$q_value
    RP.custom.out$CD.fc <- Cuffdiff.out$value_2/Cuffdiff.out$value_1
    RP.custom.out$CL.mean.FPKM1 <- rowMeans(RP.expr.data(cell.line.in, treatment1.in, treatment2.in, day.in)[,1:3])
    RP.custom.out$CL.mean.FPKM2 <- rowMeans(RP.expr.data(cell.line.in, treatment1.in, treatment2.in, day.in)[,4:6])
    RP.custom.out$CD.FPKM1 <- Cuffdiff.out$value_1
    RP.custom.out$CD.FPKM2 <- Cuffdiff.out$value_2
    
    #filter - remove genes that are not expressed in at least one RP or CD output condition - threshold as func argument
    RP.custom.out <- subset(RP.custom.out, (RP.custom.out$CL.mean.FPKM1>FPKM.filtercuttoff | RP.custom.out$CL.mean.FPKM1>FPKM.filtercuttoff) & (RP.custom.out$CD.FPKM1 >FPKM.filtercuttoff | RP.custom.out$CD.FPKM1 >FPKM.filtercuttoff) )
    print("genes remaining after FPKM filter")
    print(length(RP.custom.out$RP.pfp))
    
    #quick check for simple correlations (pearson)
    print("FPKM1 between CL and CD")
    print(cor.test(RP.custom.out$CL.mean.FPKM1,RP.custom.out$CD.FPKM1, na.rm = TRUE))
    print("FPKM2 between CL and CD")
    print(cor.test(RP.custom.out$CL.mean.FPKM2,RP.custom.out$CD.FPKM2, na.rm = TRUE))
    print("FC between RP and CD")
    print(cor.test(RP.custom.out$RP.fc,RP.custom.out$CD.fc, na.rm = TRUE))
    
    #apply FC cuttoff for each dataset
    RP.fc.RPcut <- subset(RP.custom.out,RP.custom.out$RP.fc > FC.filtercuttoff | RP.custom.out$RP.fc < 1/FC.filtercuttoff)
    print("genes remaining after fc filter on RP")
    print(length(RP.fc.RPcut$RP.pfp))
    RP.fc.CDcut <- subset(RP.custom.out,RP.custom.out$CD.fc > FC.filtercuttoff | RP.custom.out$CD.fc < 1/FC.filtercuttoff)
    print("genes remaining after fc filter on CD")
    print(length(RP.fc.CDcut$RP.pfp))
    RP.fc.Bothcut <- RP.custom.out[intersect(row.names(RP.fc.RPcut),row.names(RP.fc.CDcut)),]
    print("genes remaining after fc filter on both")
    print(length(RP.fc.Bothcut$RP.pfp))
    
    print("Number of sigdiff genes for RP @ pfp<0.05")
    print(sum(RP.fc.RPcut$RP.pfp > 0 & RP.fc.RPcut$RP.pfp < 0.05))
    print("Number of sigdiff genes for CD @ qval<0.05")
    print(sum(RP.fc.CDcut$CD.qval > 0 & RP.fc.CDcut$CD.qval < 0.05))
    print("Number sigdiff genes shared between RP and CD")
    print(sum(RP.fc.Bothcut$RP.pfp > 0 & RP.fc.Bothcut$RP.pfp < 0.05 & RP.fc.Bothcut$CD.qval > 0 & RP.fc.Bothcut$CD.qval < 0.05))
    
    return(RP.custom.out[order(RP.custom.out$RP.pfp),])
  }

#generate data.frames for all diffexpr comparisons
  MA1vME1.2 <- RP.custom("M","A","E","1")
  MA4vME4.2 <- RP.custom("M","A","E","4")
  MA5vME5.2 <- RP.custom("M","A","E","5")
  MA5vMB5.2 <- RP.custom("M","A","B","5")
  MA5vMG5.2 <- RP.custom("M","A","G","5")
  FA1vFE1.2 <- RP.custom("F","A","E","1")
  FA4vFE4.2 <- RP.custom("F","A","E","4")
  FA5vFE5.2 <- RP.custom("F","A","E","5")
  FA5vFB5.2 <- RP.custom("F","A","B","5")
  FA5vFG5.2 <- RP.custom("F","A","G","5")
#Notes: CD consistantly resulted in nearly an order of magnitude more sigdiff genes
  #correlations were good for most tests, a couple conditions had poor correlation at the FC level
  #aparent poor overlap between genes identified by the two tests
  # RP presumes "only a minority of all the features measured are upregulated (or downregulated)"
    #perhaps this holds less true for datasets with large scale (epigenetic inhibitor combinations) diffexpr changes?
  
#export outputs to .csv
  output.names <- c("MA1vME1.2","MA4vME4.2","MA5vME5.2","MA5vMB5.2","MA5vMG5.2","FA1vFE1.2","FA4vFE4.2","FA5vFE5.2","FA5vFB5.2","FA5vFG5.2")
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data Rankprod output/")
  for (x in output.names) {
    output.write <- get(x)
    write.csv(output.write, paste(x, "Rankprod output.csv"), quote = FALSE)
  } 
  rm(x)
  
  
############################################################
#for parsing in/out provided cuffdiff data
  cuffdiff.names <- c("MA1vME1","MA4vME4","MA5vME5","MA5vMB5","MA5vMG5","FA1vFE1","FA4vFE4","FA5vFE5","FA5vFB5","FA5vFG5")
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data cuffdiff output/")
  for (x in cuffdiff.names) {
    output.write <- get(x)
    write.csv(output.write, paste(x, ".csv",sep = ""), quote = FALSE)
  } 
  rm(x)
  
  #parse in from provided files
  for (x in list.files()) {
    in.name <- strsplit(x,".csv")[[1]]
    in.data <- read.csv(x,header = TRUE, stringsAsFactors = FALSE)
    assign(in.name, in.data)
  }
  rm(x,in.data,in.name)
  
  
##################################################################################################
###Network analyis################################################################################
##################################################################################################
  
#redo RankProduct analysis, but don't filter against old cuffdiff output data
RP.custom.2 <- function(cell.line.in, treatment1.in, treatment2.in, day.in,FPKM.filtercuttoff = 0,FC.filtercuttoff = 1.5){
  RP.custom.out <- data.frame(row.names = RQD_filtered.geneIDs)
  RP.custom.in <- RankProducts(RP.expr.data(cell.line.in, treatment1.in, treatment2.in, day.in),cl,logged = FALSE,gene.names = RQD_filtered.geneIDs) #random seed not fixed as above
  #print(topGene(RP.custom.in,cutoff = 0.05, method = "pfp",logged = FALSE, gene.names = RQD_filtered.geneIDs))
  print("pfp<0.05")
  print(sum(RP.custom.in$pfp <0.05))
  #data =.frame key for filtering fc values based on upregulated or downregulated
  RP.updownkey <- data.frame ("up" = (1/RP.custom.in$AveFC)>1,"down" = (1/RP.custom.in$AveFC)<1)
  RP.custom.out$RP.pfp <- rowSums(RP.updownkey*RP.custom.in$pfp)
  RP.custom.out$RP.fc <- 1/RP.custom.in$AveFC # invert so that >1 = upreg and <1 = downreg (class2/class1) as in cuffdiff output
  #fold change for IPA must be in -inf to -1 and 1 to inf format - must convert downregulated genes
  RP.custom.out$RP.fc[RP.updownkey$down,] <- -1*(1/RP.custom.out$RP.fc[RP.updownkey$down,])
  RP.custom.out$CL.mean.FPKM1 <- rowMeans(RP.expr.data(cell.line.in, treatment1.in, treatment2.in, day.in)[,1:3])
  RP.custom.out$CL.mean.FPKM2 <- rowMeans(RP.expr.data(cell.line.in, treatment1.in, treatment2.in, day.in)[,4:6])
  
  #filtering
  
  print("filter for FPKM<1 in both conditions - remaining genes:")
  RP.custom.out <- RP.custom.out[!(RP.custom.out$CL.mean.FPKM1 <1 & RP.custom.out$CL.mean.FPKM2 <1 ),]
  print(length(RP.custom.out$RP.pfp))
  print("pfp<0.05")
  print(sum(RP.custom.in$pfp <0.05))
  
  print("filter for |fc|>=1.5 in both conditions - remaining genes:")
  RP.custom.out <- RP.custom.out[abs(RP.custom.out$RP.fc) >=1.5,]
  print(length(RP.custom.out$RP.pfp))
  print("pfp<0.05")
  print(sum(RP.custom.in$pfp <0.05))
  return(RP.custom.out)
  
  }
  
  MA1vME1.2.ipa <- RP.custom.2("M","A","E","1")
  MA4vME4.2.ipa <- RP.custom.2("M","A","E","4")
  MA5vME5.2.ipa <- RP.custom.2("M","A","E","5")
  MA5vMB5.2.ipa <- RP.custom.2("M","A","B","5")
  MA5vMG5.2.ipa <- RP.custom.2("M","A","G","5")
  FA1vFE1.2.ipa <- RP.custom.2("F","A","E","1")
  FA4vFE4.2.ipa <- RP.custom.2("F","A","E","4")
  FA5vFE5.2.ipa <- RP.custom.2("F","A","E","5")
  FA5vFB5.2.ipa <- RP.custom.2("F","A","B","5")
  FA5vFG5.2.ipa <- RP.custom.2("F","A","G","5")  
  
IPAnames.2 <-  c("MA1vME1.2.ipa","MA4vME4.2.ipa","MA5vME5.2.ipa","MA5vMB5.2.ipa","MA5vMG5.2.ipa","FA1vFE1.2.ipa","FA4vFE4.2.ipa","FA5vFE5.2.ipa","FA5vFB5.2.ipa","FA5vFG5.2.ipa")
  
#Parse out data for IPA network analysis
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data for IPA import/")  
  
  for (x in IPAnames.2) {
    data <- get(x)
    print(x)
    write.table(data,file = paste(x,".txt"),sep = "\t",quote = FALSE,row.names = TRUE,col.names = FALSE)
  }
  rm(x,data)
  
#data imported to IPA with pfp<0.05 as cutoff
#insuficient number of genes for analysis of day 1 data for either cell line
#data manually exported to tab deliited text files and are parsed in below
#no GO analysis for the purposes of this
  
  #read in tab delim files
  #pathway
  #set wd to pathway analysis files
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data for IPA output/Pathway/")  
  for(x in list.files()){
    xx <- unlist(strsplit(x,split = ".txt"))
    y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2)
    #get molecules as list instead of comma delim
    y$list <- sapply(y$Molecules, function(x){unlist(strsplit(x,split = ","))})
    y <- y[,c("Ingenuity.Canonical.Pathways","X.log.p.value.","z.score","list")]
    assign(paste(xx,"IPA_pathway","3",sep = "_"),y)
  }
  rm(x,xx,y)
  
  #upstream
  #set wd to upstream analysis files
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data for IPA output/UR/")  
  for(x in list.files()){
    print(x)
    xx <- unlist(strsplit(x,split = ".txt"))
    
    y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2,fill = TRUE) 
    #must manually add null activation Z score cols for both A5vB5's (column does not exist in output file)
    if (xx == "MA5vMB5" | xx == "FA5vFB5" ){
      y$Activation.z.score <- NA
    }
    #get molecules as list instead of comma delim
    y$list <- sapply(y$Target.molecules.in.dataset, function(x){unlist(strsplit(x,split = ","))})
    y <- y[,c("Upstream.Regulator","p.value.of.overlap","Activation.z.score","list")]
    assign(paste(xx,"IPA_upstream","3",sep = "_"),y)
  }
  rm(x,xx,y)
  
  #make named list of df names for each
  IPA.all <- c("FA4vFE4", "FA5vFB5", "FA5vFE5", "FA5vFG5", "MA4vME4", "MA5vMB5", "MA5vME5", "MA5vMG5")
  IPA.pathway.3 <- sapply(IPA.all,function(x) {paste(x,"IPA_pathway_3",sep = "_")})
  IPA.upstream.3 <- sapply(IPA.all,function(x) {paste(x,"IPA_upstream_3",sep = "_")})
  
########### visualization via pheatmap ##################################################
  
  #pheatmap
  install.packages("pheatmap") #installs rcolorbrewer
  library("pheatmap")

  
#new pheatmap generating function that compares IPA outputs from old and RP data
    
  #accepts list (char) of regulators or pathways in order - plots heat map forthose
  #use two color scales (trick into displaying both by hybrid color scale - add 100 to z-scores)
  #2 types accepted: "IPA_pathway" ; "IPA_upstream2"
  # p-val is column 2 for all three; z-score is column 3 for all three
  #grey indicates no score was provided to distinguish from a non-significant score (white)
  #differencs from IPA.heatmap3
    #type is only pathway or UR, internally aquires both version 2 (old data) and 3 (Rankproduct data)
    #dont incorperate UTC genes
  IPA.heatmap4 <- function(IPA.id.list,type, width = 15, height = 15, analysis.number = "3"){ 
    
    DF<-data.frame(row.names = IPA.id.list)#set empty DF to be filled with heatmap data
    
    names<-  c(sapply(IPA.all,function(x) {paste(x,type,"2",sep = "_")}),sapply(IPA.all,function(x) {paste(x,type,analysis.number,sep = "_")}))
    #re-order to place cell lines next to each other
    names<- names[c(1:4,9:12,5:8,13:16)]
    for(x in names){
      DF.IPA <- get(x)
      DF.IPA[is.na(DF.IPA)] <- 50 #!!!!!!!!changed from origional - change NaN to 50 not 0 (to distinguish on color scale)
      DF.ipa.test.2 <<-DF.IPA #for troubleshooting
      if(type == "IPA_upstream"){ #subsets for p<0.05
        DF.IPA <- DF.IPA[(DF.IPA[,2] < 0.05),] #exlude non-sig (exclude p>0.05)
      }else{
        DF.IPA <- DF.IPA[(DF.IPA[,2] > (-log10(0.05))),] #exlude non-sig (exclude p>0.05)
      }
      
      row.names(DF.IPA)<- DF.IPA[,1]
      
      if(type == "IPA_upstream"){ #converts to -log p_val for upstream to match pathway
        DF[,paste(x,"-log_p_value",sep = " ")] <- -log10(DF.IPA[IPA.id.list,2])
      }else{
        DF[,paste(x,"-log_p_value",sep = " ")] <- DF.IPA[IPA.id.list,2]
      }
      
      DF[,paste(x,"activation_z_score",sep = " ")]<-DF.IPA[IPA.id.list,3] + 100
      DF[is.na(DF)]<-0 # all NA's to 0'
    }
    
    DF.ipa.test <<- DF #assign as global variable for troubleshooting
    
    #set color parameters
    hm.breaks <- c(0,seq(-log10(0.05),max(DF[DF<75]),length = 10)) #sets breaks for color scale, total # of breaks = 31 
    #if no negative z scores in entire data frame exclude those breaks/colors
    if(length(DF[DF>75 & DF<100]) != 0){ 
      hm.breaks <- c(hm.breaks,seq(min(DF[DF>75])-0.00001,max(DF[DF<98]),length = 10))
    }else{
      hm.breaks <- c(hm.breaks,99.9999)
    }
    hm.breaks <- c(hm.breaks,seq(min(DF[DF>102]),max(DF[(DF>102 & DF<150)]),length = 10),151)
    print("color scale breaks:")
    print(hm.breaks)
    
    col1 <- (colorRampPalette(c("black",'red'))(9+4))[5:13]
    col2 <- (colorRampPalette(c("purple",'white'))(9+4))[1:9]
    col3 <- (colorRampPalette(c("white",'orange'))(9+4))[5:13]
    hm.colors <- c("white",col1,"grey") #should be 1 less than hm.breaks; white = no enrichment
    #if no negative z scores in entire data frame exclude those breaks/colors
    if(length(DF[DF>75 & DF<100]) != 0){
      hm.colors <- c(hm.colors,col2)
    }
    hm.colors <- c(hm.colors,"white",col3,"grey")
    
    #make annotation data frame
    hm.col.annotation <- data.frame(Cell.Line <- factor(rep(c("FLAM76","MMM1"),each = 16)),
                                    Day <- factor(rep(c("Day 4","Day 4",rep("Day 5.5",6)),4)),
                                    EPZ <- factor(rep(c("+","+","+","+","-","-","+","+"),4)),
                                    PAN <- factor(rep(c("-","-","-","-","+","+","+","+"),4)),
                                    METHOD <- factor(rep(c("cuffdiff","rankproduct"),each =8, times =2)),
                                    row.names = colnames(DF))
    colnames(hm.col.annotation)<-c("Cell Line","Sampling Day     ","EPZ6438","Panobinostat","Method")
    #specify colors for annotation
    hm.col.annotation.colors <- list("Cell Line" = c(FLAM76 = "gold",MMM1 ="darkorange2"),
                                     "Sampling Day     " = c("Day 4" = "royalblue","Day 5.5" = "navy"),
                                     "EPZ6438" = c("+" = "black","-" = "white"),
                                     "Panobinostat" = c("+" = "black","-" = "white"),
                                     "Method" = c("cuffdiff" = "red","rankproduct" = "green"))
    
    #generate heatmap
    if(height < 1) {hm.ylab <- FALSE}else{hm.ylab <- TRUE} #dont display genes if cells are too short
    pheatmap(DF,
             color = hm.colors,
             breaks = hm.breaks,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             gaps_col = c(2,2,2,4,6,8,8,8,8,8,8,10,10,10,12,14,16,16,16,16,16,16,18,18,20,22,24,24,24,24,24,24,26,26,26,28,30),#by day
             cellwidth = width,
             cellheight = height,
             treeheight_row = 0,
             annotation_col = hm.col.annotation,
             annotation_colors = hm.col.annotation.colors,
             show_colnames = FALSE,
             show_rownames = hm.ylab,
             legend_labels = c("log2 Fold change"),
             fontsize_row = height +3
    )
  }
  
  IPA.final.upstream.list <-c(
    "EZH2",
    "Hdac",
    "HDAC1",
    "HDAC2",
    "CCND1",
    "TP53",
    "MYC",
    "MYCN",
    "BTK",
    "MAPK1",
    "ERBB2",
    "TNF",
    "ATM",
    "Jnk",
    "BRD4",
    "P38 MAPK",
    "EIF2AK2",
    "SMARCA4",
    "SMARCB1",
    "Interferon alpha",
    "IFNG",
    "STAT1",
    "IRF4",
    "SOCS1",
    "SOCS3",
    "estrogen receptor",
    "ESR1",
    "BCL6",
    "PTEN",
    "HIF1A",
    "TGFB1",
    "NFkB (complex)",
    "JMJD6"
  )
  
  
  IPA.final.pathway.list <-c(
    "Interferon Signaling",
    "B Cell Development",
    "Antigen Presentation Pathway",
    "Molecular Mechanisms of Cancer",
    "Calcium-induced T Lymphocyte Apoptosis",
    "IL-4 Signaling",
    "Th1 and Th2 Activation Pathway",
    "Integrin Signaling",
    "Phagosome Formation",
    "Axonal Guidance Signaling",
    "Epithelial Adherens Junction Signaling",
    "IL-8 Signaling",
    "Mitochondrial Dysfunction",
    "Endoplasmic Reticulum Stress Pathway",
    "Granzyme A Signaling",
    "Estrogen-mediated S-phase Entry",
    "Cyclins and Cell Cycle Regulation",
    "Cell Cycle Control of Chromosomal Replication",
    "Fatty Acid β-oxidation I",
    "Sirtuin Signaling Pathway",
    "tRNA Charging"
  )
  
  #generate heatmaps
  IPA.heatmap4(IPA.final.upstream.list,type = "IPA_upstream",height = 15,width = 15)
  IPA.heatmap4(IPA.final.pathway.list,type = "IPA_pathway",height = 15,width = 15)
  
  
  
#compile new lists of regulators/pathways based on most frequent/significant RP hits
  #test that RP hits are also contained in cuffdiff data
  
  #iterate through each of the rankproduct anallysis ipa output files and generate a list of most frequently represented hits
  
    #upstream regulators
    IPA.upstream.3.unique <- character()
    for (x in IPA.upstream.3) {
      IPA.upstream.3.unique <- unique(c(IPA.upstream.3.unique,get(x)$Upstream.Regulator))
    }
    rm(x)
    
    IPA.upstream.3.freqcount <- data.frame(row.names = IPA.upstream.3.unique)
    for (x in IPA.upstream.3) {
      print(x)
      data <- get(x)
      index <- rownames(IPA.upstream.3.freqcount)[which(rownames(IPA.upstream.3.freqcount) %in% data$Upstream.Regulator)]
      print(length(index))
      row.names(data)<-data$Upstream.Regulator
      IPA.upstream.3.freqcount[index,x] <- data[index,"p.value.of.overlap"]
    }
    rm(x,data,index)
    #remove insignificant values before freq count (IPA does kick back some)
    IPA.upstream.3.freqcount[IPA.upstream.3.freqcount >0.05]<- NA
    IPA.upstream.3.freqcount$freqcount <- rowSums(!(is.na(IPA.upstream.3.freqcount)))
    #sort by freq# and then lowest mena pval
    IPA.upstream.3.freqcount <- IPA.upstream.3.freqcount[order(-IPA.upstream.3.freqcount$freqcount, rowMeans(IPA.upstream.3.freqcount[,1:8], na.rm = TRUE)),]
    
    #regulators for heatmap
    IPA.RP.upstream.list <- rownames(IPA.upstream.3.freqcount)[1:100]
  #visualize heatmap
  IPA.heatmap4(IPA.RP.upstream.list,type = "IPA_upstream",height = 6,width = 6)
    
    
    
    
    
    
    #canonical pathways
    IPA.pathway.3.unique <- character()
    for (x in IPA.pathway.3) {
      IPA.pathway.3.unique <- unique(c(IPA.pathway.3.unique,get(x)$Ingenuity.Canonical.Pathways))
    }
    rm(x)
    
    IPA.pathway.3.freqcount <- data.frame(row.names = IPA.pathway.3.unique)
    for (x in IPA.pathway.3) {
      print(x)
      data <- get(x)
      index <- rownames(IPA.pathway.3.freqcount)[which(rownames(IPA.pathway.3.freqcount) %in% data$Ingenuity.Canonical.Pathways)]
      print(length(index))
      row.names(data)<-data$Ingenuity.Canonical.Pathways
      IPA.pathway.3.freqcount[index,x] <- data[index,"X.log.p.value."]
    }
    rm(x,data,index)
    #remove insignificant values before freq count (IPA does kick back some)
    IPA.pathway.3.freqcount[IPA.pathway.4.freqcount < -log10(0.05)]<- NA
    IPA.pathway.3.freqcount$freqcount <- rowSums(!(is.na(IPA.pathway.3.freqcount)))
    #sort by freq# and then lowest mena pval
    IPA.pathway.3.freqcount <- IPA.pathway.3.freqcount[order(-IPA.pathway.3.freqcount$freqcount, -rowMeans(IPA.pathway.3.freqcount[,1:8], na.rm = TRUE)),]
    
    #regulators for heatmap
    IPA.RP.pathway.list <- rownames(IPA.pathway.3.freqcount)[1:100]
  #visualize heatmap
  IPA.heatmap4(IPA.RP.pathway.list,type = "IPA_pathway",height = 6,width = 6)

  
############ additional IPA analysis from RP data #############################################################
  #open question - the lack of overlap between informatic results is expected to a point do to 
    #the differences in number of sigdiff genes output by each method
    #are these diferences reflecting a dissagreement on which genes are differentially expressed 
    #due to the method only producing a smaller number of differentially expressed genes
    #...or are the distributions in significance correlate such that RP method simply produces a much
    #more conservative threshold?
  #test: lower cuttoff from pf<0.05 to pfp<0.15 (roughly doubles genes considered by IPA as significant)
      #~600 in combos, ~300 in single agents at day 5 (except FA5vFB5 where very little change is seen in any tested method)
    #use same fc cutoffs (this was filtered before IPA inport anyway, pfp cuttoff was not)
    #exclude day one observation
  
  
  #parse in new IPA output files (tab delim files manually downloaded in IPA client)
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data for IPA output 2/Pathway/")  
  for(x in list.files()){
    xx <- unlist(strsplit(x,split = ".txt"))
    y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2)
    #get molecules as list instead of comma delim
    y$list <- sapply(y$Molecules, function(x){unlist(strsplit(x,split = ","))})
    y <- y[,c("Ingenuity.Canonical.Pathways","X.log.p.value.","z.score","list")]
    assign(paste(xx,"IPA_pathway","4",sep = "_"),y)
  }
  rm(x,xx,y)
  
  #upstream
  #set wd to upstream analysis files
  setwd("D:/Post-doc applications/TG at UCSF/Re-quantfied EZH2 paper RNASEQ data for IPA output 2/UR/")  
  for(x in list.files()){
    print(x)
    xx <- unlist(strsplit(x,split = ".txt"))
    
    y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2,fill = TRUE) 
    #must manually add null activation Z score cols for both A5vB5's (column does not exist in output file)
    if (xx == "MA5vMB5" | xx == "FA5vFB5" ){
      y$Activation.z.score <- NA
    }
    #get molecules as list instead of comma delim
    y$list <- sapply(y$Target.molecules.in.dataset, function(x){unlist(strsplit(x,split = ","))})
    y <- y[,c("Upstream.Regulator","p.value.of.overlap","Activation.z.score","list")]
    assign(paste(xx,"IPA_upstream","4",sep = "_"),y)
  }
  rm(x,xx,y)
  
  #make named list of df names for each
  IPA.all <- c("FA4vFE4", "FA5vFB5", "FA5vFE5", "FA5vFG5", "MA4vME4", "MA5vMB5", "MA5vME5", "MA5vMG5")
  IPA.pathway.4 <- sapply(IPA.all,function(x) {paste(x,"IPA_pathway_4",sep = "_")})
  IPA.upstream.4 <- sapply(IPA.all,function(x) {paste(x,"IPA_upstream_4",sep = "_")})
  
  
  #generate heatmaps
  #added analysis.number argument to heatmap func (char: expects "3" or "4" where 4 is ipa output with more liberal cuttoff)
    #default = "3" for back compatability when running above code
  IPA.heatmap4(IPA.final.upstream.list,type = "IPA_upstream",height = 15,width = 15, analysis.number = "4")
  IPA.heatmap4(IPA.final.pathway.list,type = "IPA_pathway",height = 15,width = 15, analysis.number = "4")
  
  #iterate through each of the rankproduct anallysis ipa output files and generate a list of most frequently represented hits
  
  #upstream regulators
  IPA.upstream.4.unique <- character()
  for (x in IPA.upstream.4) {
    IPA.upstream.4.unique <- unique(c(IPA.upstream.4.unique,get(x)$Upstream.Regulator))
  }
  rm(x)
  
  IPA.upstream.4.freqcount <- data.frame(row.names = IPA.upstream.4.unique)
  for (x in IPA.upstream.4) {
    print(x)
    data <- get(x)
    index <- rownames(IPA.upstream.4.freqcount)[which(rownames(IPA.upstream.4.freqcount) %in% data$Upstream.Regulator)]
    print(length(index))
    row.names(data)<-data$Upstream.Regulator
    IPA.upstream.4.freqcount[index,x] <- data[index,"p.value.of.overlap"]
  }
  rm(x,data,index)
  #remove insignificant values before freq count (IPA does kick back some)
  IPA.upstream.4.freqcount[IPA.upstream.4.freqcount >0.05]<- NA
  IPA.upstream.4.freqcount$freqcount <- rowSums(!(is.na(IPA.upstream.4.freqcount)))
  #sort by freq# and then lowest mena pval
  IPA.upstream.4.freqcount <- IPA.upstream.4.freqcount[order(-IPA.upstream.4.freqcount$freqcount, rowMeans(IPA.upstream.4.freqcount[,1:8], na.rm = TRUE)),]
  
  #regulators for heatmap
  IPA.RP.4.upstream.list <- rownames(IPA.upstream.4.freqcount)[1:100]
  #visualize heatmap
  IPA.heatmap4(IPA.RP.4.upstream.list,type = "IPA_upstream",height = 6,width = 6, analysis.number = "4")
  
  
  
  
  
  
  #canonical pathways
  IPA.pathway.4.unique <- character()
  for (x in IPA.pathway.4) {
    IPA.pathway.4.unique <- unique(c(IPA.pathway.4.unique,get(x)$Ingenuity.Canonical.Pathways))
  }
  rm(x)
  
  IPA.pathway.4.freqcount <- data.frame(row.names = IPA.pathway.4.unique)
  for (x in IPA.pathway.4) {
    print(x)
    data <- get(x)
    index <- rownames(IPA.pathway.4.freqcount)[which(rownames(IPA.pathway.4.freqcount) %in% data$Ingenuity.Canonical.Pathways)]
    print(length(index))
    row.names(data)<-data$Ingenuity.Canonical.Pathways
    IPA.pathway.4.freqcount[index,x] <- data[index,"X.log.p.value."]
  }
  rm(x,data,index)
  #remove insignificant values before freq count (IPA does kick back some)
  IPA.pathway.4.freqcount[IPA.pathway.4.freqcount < -log10(0.05)]<- NA
  IPA.pathway.4.freqcount$freqcount <- rowSums(!(is.na(IPA.pathway.4.freqcount)))
  #sort by freq# and then lowest mena pval
  IPA.pathway.4.freqcount <- IPA.pathway.4.freqcount[order(-IPA.pathway.4.freqcount$freqcount, -rowMeans(IPA.pathway.4.freqcount[,1:8], na.rm = TRUE)),]
  
  #regulators for heatmap
  IPA.RP.4.pathway.list <- rownames(IPA.pathway.4.freqcount)[1:100]
  #visualize heatmap
  IPA.heatmap4(IPA.RP.4.pathway.list,type = "IPA_pathway",height = 6,width = 6, analysis.number = "4")
  
  
  #generate table to illustrate overlapping gene sets for each condition
    #two tables, one for RP data at two threasholds
  
  final.row.names <- c(
    "FLAM76 DAY 4 EPZ6438",
    "FLAM76 DAY 5 Panobinostat",
    "FLAM76 DAY 5 EPZ6438",
    "FLAM76 DAY 5 Combination",
    "MMM1 DAY 4 EPZ6438",
    "MMM1 DAY 5 Panobinostat",
    "MMM1 DAY 5 EPZ6438",
    "MMM1 DAY 5 Combination"
  )
  
  
  #declare data.frame to populate with comparisons
    #IPA.all = "FA4vFE4" "FA5vFB5" "FA5vFE5" "FA5vFG5" "MA4vME4" "MA5vMB5" "MA5vME5" "MA5vMG5"
  #RP pfp<0.05
  CD.RP.compare.1 <- data.frame(row.names = IPA.all)
  for (x in IPA.all) {
    data <- get(x)
    #filter for |fc|>1.5, FPKM>1 in at least one (treatment/control) and qval<0.05
    data <- data[abs(data$log2.fold_change.)>log2(1.5) & (data$value_1 >=1 | data$value_2 >= 1) & data$q_value < 0.05,]
    CD.up <- data$gene[data$log2.fold_change.>0]
    CD.down <- data$gene[data$log2.fold_change.<0]
    CD.RP.compare.1[x,"Cuffdiff Upregulated"] <- length(CD.up)
    CD.RP.compare.1[x,"Cuffdiff Downregulated"] <- length(CD.down)
    
    data <- get(paste(x,".2.ipa", sep = ""))
    #already filtered for FPKM and FC (same cutoffs as CD data)
    #filter for pfp cutoff @ 0.05
    data <- data[data$RP.pfp <0.05,]
    RP.up <- row.names(data)[data$RP.fc>0]
    RP.down <- row.names(data)[data$RP.fc<0]
    CD.RP.compare.1[x,"Rankproduct Upregulated"] <- length(RP.up)
    CD.RP.compare.1[x,"Rankproduct Downregulated"] <- length(RP.down)
    
    CD.RP.compare.1[x,"Overlap Upregulated"] <- length(intersect(CD.up,RP.up))
    CD.RP.compare.1[x,"Overlap Downregulated"] <- length(intersect(CD.down,RP.down))
  }
  rm(x,data)
  row.names(CD.RP.compare.1)<- final.row.names
  
  
  #RP pfp<0.15
  CD.RP.compare.2 <- data.frame(row.names = IPA.all)
  for (x in IPA.all) {
    data <- get(x)
    #filter for |fc|>1.5, FPKM>1 in at least one (treatment/control) and qval<0.05
    data <- data[abs(data$log2.fold_change.)>log2(1.5) & (data$value_1 >=1 | data$value_2 >= 1) & data$q_value < 0.05,]
    CD.up <- data$gene[data$log2.fold_change.>0]
    CD.down <- data$gene[data$log2.fold_change.<0]
    CD.RP.compare.2[x,"Cuffdiff Upregulated"] <- length(CD.up)
    CD.RP.compare.2[x,"Cuffdiff Downregulated"] <- length(CD.down)
    
    data <- get(paste(x,".2.ipa", sep = ""))
    #already filtered for FPKM and FC (same cutoffs as CD data)
    #filter for pfp cutoff @ 0.05
    data <- data[data$RP.pfp <0.15,]
    RP.up <- row.names(data)[data$RP.fc>0]
    RP.down <- row.names(data)[data$RP.fc<0]
    CD.RP.compare.2[x,"Rankproduct Upregulated"] <- length(RP.up)
    CD.RP.compare.2[x,"Rankproduct Downregulated"] <- length(RP.down)
    
    CD.RP.compare.2[x,"Overlap Upregulated"] <- length(intersect(CD.up,RP.up))
    CD.RP.compare.2[x,"Overlap Downregulated"] <- length(intersect(CD.down,RP.down))
  }
  rm(x,data)
  row.names(CD.RP.compare.2)<- final.row.names
  
  #save comparison data.frames to .csv files
  setwd("D:/Post-doc applications/TG at UCSF/")
  write.csv(CD.RP.compare.1, "Overlapping genes between methods.csv")
  write.csv(CD.RP.compare.2, "Overlapping genes between methods loosened cutoff.csv")
  