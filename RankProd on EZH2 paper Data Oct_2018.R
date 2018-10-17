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
  
  
############################################################333
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