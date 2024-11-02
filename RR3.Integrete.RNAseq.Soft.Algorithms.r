
####### 
## source("/home/kes21/_Scripts/00..MAIN.MAIN.DIRs.r");  source(paste(MAIN.MAIN.DIR...Pipeline..RNseq,"_ghub_Talal/RR3.Integrete.RNAseq.Soft.Algorithms.r",sep=""))

#######  
## /n/groups/markianos/klaus/88_Programs/R_Versions/R-3.6.3/bin/R

####################################################################################################       
print("##################################################################")
print("########          Integrete RNAseq Soft Algorithms        ########") 
print("#                                                                #")
print("#   Klaus Schmitz Abe (Klaus.Schmitz-Abe@childrens.harvard.edu)  #")
print("##################################################################")
rm(list=ls());	print(date());   print("");      Initial_TIME= date();   WD.ORG= getwd(); 

####################################################################################################   
####################################################################################################    
############ missing

####################################################################################################  
####################################################################################################  
Using_argumnets_in_R = 1
if (Using_argumnets_in_R==1) {  Inpu_args = commandArgs(trailingOnly=T)
	########################
	Project.RNA..run                 = as.character(Inpu_args[ 1]);
	DIR_Annotations.Source.hg..VEC.TX= as.character(Inpu_args[ 2]);
	print..RNAseq_results			 = as.integer  (Inpu_args[ 3]);	   
	print..RNAseq..script            = as.integer  (Inpu_args[ 4]); 
	######
	Run..Lima.models..print_tests    = as.integer  (Inpu_args[ 5]);
	Remove.DIR...Results.Models.Step3= as.integer  (Inpu_args[ 6]);
	RNAseq..Project_num.Ver          = as.character(Inpu_args[ 7]);
	RUN..Trimmomatic..ONE  		     = as.character(Inpu_args[ 8]);
	########################
	DIR_Annotations.Source.hg..VEC= unlist(strsplit(DIR_Annotations.Source.hg..VEC.TX, ";;;", fixed=TRUE))
	########################
	Inputs.Bsub= 1;			Runnning_testing_functions_inputs= c(0,1,2)[1];
	print(paste("INPUTS --> ",paste(Inpu_args,collapse=", ")));  
}
if( is.na(Project.RNA..run[1]) ) {
	############
	Inputs.Bsub= 0;			Runnning_testing_functions_inputs= c(0,1,2)[1];    #### 1: stop at the end;   2= run testing_functions
	######
    Project.RNA..run= c("Example..P1")[1]
    ######
    DIR_Annotations.Source.hg..VEC= c("Ensembl:Ch37","NCBI:Ch37","UCSC:hg19","UCSC:hg38","Ensembl:Ch38","RSEM","Others")[5]; 		RUN..Trimmomatic..ONE = c("No","Yes")[2] 
	############ 
	RNAseq..Project_num.Ver= c("AA","BB","CC")[1];		Remove.DIR...Results.Models.Step3= c(1,0)[1]; 
	######
	print..RNAseq_results=c(0,1,2)[1];		   print..RNAseq..script=c(0,1,2)[1];		Run..Lima.models..print_tests= c(1,0)[1]   
}
####################################################################################################################
####################################################################################################################
MAIN_WROKING_Directory= paste("/home/kes21/.configset/.emacs.dir2/XXX_Loco2/00_Pipeline..VExP/RR0.RNAseq/_ghub_Talal/",sep="")

####################################################################################################################
####################################################################################################################
CASE..Res.Alg= c(1,2,3,4)[2];   	CASE..Res.GO= c(1,2,3,4)[2];          				RNAseq..FORCE_USE.only.one.libType..VEC= c(0,1)[1];	  	
################
MIN.logFC..Filters = c( -1, 1.2, 1.5, 2)[CASE..Res.Alg];		MAX.PValue..Filters= c(-1, 0.05, 0.05, 0.01)[CASE..Res.Alg];		MAX.FDR..Filters1= c(-1, 0.05, 0.05, 0.01)[CASE..Res.Alg];		MAX.PValue.Path..Filters1= c(-1, 1e-04, 1e-05, 1e-07)[CASE..Res.Alg]; 						         
MIN.logFC..Res.GO  = c( -1, 1.2, 1.5, 2)[CASE..Res.GO ]; 		MAX.PValue..Res.GO = c(-1, 0.05, 0.05, 0.01)[CASE..Res.GO ];		MAX.FDR..Filters2= c(-1,  0.1, 0.1 , 0.05)[CASE..Res.Alg];		MAX.PValue.Path..Filters2= c(-1,  0.05, 1e-03, 1e-04)[CASE..Res.Alg];   
COL.MAX.MEAN..Res.GO= c(".Max",".Mean")[1];    	                                                                                                                                		   	
#########
MAX.FDR..Res.GO = c(-1, 0.05,0.05, 0.01)[CASE..Res.GO ];		MAX.FDR..Res.GO..Plot.QC = MAX.FDR..Res.GO;
#########
CASEs.PLOT..PAthway.GO.KEGG = c(1,2,3,4)[c(2:2)];				CASEs.PLOT..Overlapping.Gene= c(1,2,3,4)[c(2:3)];					PLOT.VennDiagram..GO.KEGG= c(1,0)[1];							PLOT.barcodeplot..GO.KEGG= c(1,0)[1];		  
#########
MAX.number.plots.barcodeplot= c(3,10,100)[3];        			MAX.PValue.Path..Barcodeplot= c(-1,1e-03,1e-04,1e-05)[3];			MAX.PValue.Path.KEGG.VenDiagram..VEC=c(-1,1e-02,1e-04,1e-05)[2:3];	   	
################
RNAseq..QC..filtering..Counts= c(1,0)[1];     					RNAseq..Models..VOOM= c(1,2,3)[1];									filtering.Boris1a.Boris1b.Boris1.Limma.none= c(1,2,3,4,0)[1];

################################
################################
RNAseq..Models..MIN.log2.fold= MIN.logFC..Filters;    RNAseq..Models..MIN.p.value= MAX.PValue..Filters;     RNAseq..Models..MIN.Fil1.FDR= MAX.FDR..Filters1;     RNAseq..Models..MIN.Fil2.FDR= MAX.FDR..Filters2
RNAseq..Figure..MIN.log2.fold= MIN.logFC..Filters;    RNAseq..Figure..MIN.p.value= MAX.PValue..Filters;     RNAseq..Figure..MIN.Fil1.FDR= MAX.FDR..Filters1;     RNAseq..Figure..MIN.Fil2.FDR= MAX.FDR..Filters2

####################################################################################################################
#################################################################################################################### test
#Remove.DIR...Results.Models.Step3= c(1,0)[2];			
#CASEs.PLOT..Overlapping.Gene= c(1,2,3,4)[c(4)];     PLOT.VennDiagram..GO.KEGG= c(1,0)[2];   
#PLOT.barcodeplot..GO.KEGG= c(1,0)[1];   			 MAX.PValue.Path..Barcodeplot = c(-1, 1e-04, 1e-05, 1e-07)[4];  MAX.number.plots.barcodeplot= c(3,10,100)[2]

####################################################################################################################
####################################################################################################################
PRINT.FDR.Lima=0.05;					PRINT.FDR.DESeq= 0.05;    				PRINT.FDR.edgeR= 0.05;					PRINT.PValue=0.05
#########
Run.Model..LIMMA= c(1,0)[1];			Run.Model..DESeq= c(1,0)[1];			Run.Model..DESeq2= c(1,0)[1];        Run.Model..edgeR= c(1,0)[1];			Run.edgeR..Standard.Boris..VEC=c(1,0);	      Run.Models..Plotting= c(1,0)[1];	
Run.Models..ORG.BORIS= c(0,1)[1];		Run.Models..Plotting.TEST= c(1,0)[2];	Saving..TABLE.Percentage= c(1,0)[2]				

####################################################################################################################
####################################################################################################################
MAIN.MAIN.DIR...Pipeline..RNseq  = MAIN_WROKING_Directory
MAIN.MAIN.DIR...INPUTS..RNseq    = MAIN_WROKING_Directory
MAIN.MAIN.DIR...OUTPUT..RNseq    = MAIN_WROKING_Directory
#####
MAIN.MAIN.DIR...Programs.ALL     = MAIN_WROKING_Directory 
MAIN.MAIN.DIR...Inp.DataBases.ALL= MAIN_WROKING_Directory
########################
Running_Step_3_Models_Process.Data=1;
if(Run..Lima.models..print_tests!=1){ 
	Remove.DIR...Results.Models.Step3= 0;   Running_Step_3_Models_Process.Data=0;    print..RNAseq_results=c(0,1,2)[1]
	DIR_Annotations.Source.hg..VEC= DIR_Annotations.Source.hg..VEC[1] 
	RNAseq..FORCE_USE.only.one.libType..VEC= RNAseq..FORCE_USE.only.one.libType..VEC[1]
}
#########################
ERROR.NOT.find.genes_in_path.ALL= c();   SAMPLES.QC.ERROR..ALL= c()
for(RNAseq..FORCE_USE.only.one.libType in RNAseq..FORCE_USE.only.one.libType..VEC) for(DIR_Annotations.Source.hg in DIR_Annotations.Source.hg..VEC) {
	print(paste("##################################################################################", date() ))
	#########
	## source("/home/kes21/_Scripts/00..MAIN.MAIN.DIRs.r");
	#########
	setwd(MAIN_WROKING_Directory)
	source(paste("RR1.RNAseq..Definitions.r",sep=""))
	source(paste("RR1.RNAseq..Tests.r"      ,sep=""))

	#################
	if(length(RNAseq..Models..NOT_RUN)>0) for(MOD.k in RNAseq..Models..NOT_RUN) {
		if(MOD.k=="Lima"  ) Run.Model..LIMMA = c(1,0)[2];
		if(MOD.k=="edgeR" ) Run.Model..edgeR = c(1,0)[2];
		if(MOD.k=="DESeq" ) Run.Model..DESeq = c(1,0)[2];
	}
	
	#################
	conta.test=0
	for(test_RUN.Model in RNAseq..test_run..VEC) { 
		################################## 
		conta.test= conta.test + 1;   if(conta.test>0) Remove.DIR...Results.Models.Step3= 0
		##################################
		TXT.Design.x="";       if(RNAseq..FORCE_USE.only.one.libType!=0) TXT.Design.x= paste("-L",RNAseq..FORCE_USE.only.one.libType,sep=""); 
		txt..Trimmomatic="";   if(RUN..Trimmomatic..ONE=="Yes") txt..Trimmomatic=".T" 
		#############
		FILE.Join.Algorithms..csv = paste(sub("-TestXXX",paste("-Test",test_RUN.Model,sep=""),TXT.filtering.VOOM..SVersion,fixed=TRUE),TXT.Design.x,txt..Trimmomatic,".csv",sep="");   if(Run..Lima.models..print_tests==1){ print(paste(" running for -->", FILE.Join.Algorithms..csv));    print(paste("             -->", DIR_RNAseq..Results.Models)) }
		##################################
		FILE.Join.Algorithms..Rdata = sub(".csv",".Rdata"                            ,FILE.Join.Algorithms..csv,fixed=TRUE)
		FILE.Join.Counts..csv       = sub(".csv","..A.Counts.csv"                    ,FILE.Join.Algorithms..csv,fixed=TRUE);
		FILE.Join.Percentage..csv   = sub(".csv","..A.Counts..Per.csv"               ,FILE.Join.Algorithms..csv,fixed=TRUE);  	FILE.Join.Percentage..csv= paste("XX..",FILE.Join.Percentage..csv,sep="")
		FILE.Join.Algorithms.Com.csv= sub(".csv","..B.Diff.Gene_expression.csv"      ,FILE.Join.Algorithms..csv,fixed=TRUE);
		FILE.Join.Algorithms.Fil.csv= sub(".csv","..B.Diff.Gene_exp..Filter.csv"     ,FILE.Join.Algorithms..csv,fixed=TRUE);
		FILE.Join.Algorithms..PDF2  = sub(".csv","..C.GO-KEGG_path_analysis.pdf"     ,FILE.Join.Algorithms..csv,fixed=TRUE);
		FILE.Join.TOGO..csv         = sub(".csv","..C.GO-KEGG_path_analysis.csv"     ,FILE.Join.Algorithms..csv,fixed=TRUE);
		FILE.Join.TOGO..csv.Per     = sub(".csv","..C.GO-KEGG_path_analysis..Per.csv",FILE.Join.Algorithms..csv,fixed=TRUE);  	FILE.Join.TOGO..csv.Per  = paste("XX..",FILE.Join.TOGO..csv.Per  ,sep="")
		FILE.Join.Algorithms..PDF1  = sub(".csv","..D.All_plots.pdf"                 ,FILE.Join.Algorithms..csv,fixed=TRUE);    FILE.PDF..TEST= paste("XXX..",sub(".pdf","..TEST.pdf",FILE.Join.Algorithms..PDF1,fixed=TRUE),sep="")
		FILE.Join.Algorithms..PDF3  = sub(".csv","..E.Plots.pdf"                     ,FILE.Join.Algorithms..csv,fixed=TRUE);
		
		##################################
		if(Run.Models..Plotting.TEST==1){ setwd(DIR_RNAseq..Results.Project);      pdf(FILE.PDF..TEST) }
		######
		ERROR.NOT.find.genes_in_path=c()

############################################################################################################
############################################################################################################ 
print(""); print("############################################################# Loading inputs")
#########
setwd(DIR_RNAseq..Results.RData);    load(FILE_RNAseq..Counts.Rdata..Step2);    if(Run..Lima.models..print_tests==1){ print(paste("  Load Counts      -->",nrow(Counts.ORG), FILE_RNAseq..Counts.Rdata..fastq )) }
setwd(DIR_RNAseq..Inputs.Project);   Targets.Full= read.csv(FILE.Samples_Info, header = TRUE, stringsAsFactors = FALSE);	 if(Run..Lima.models..print_tests==1){ print(paste("  Load samples id  -->",FILE.Samples_Info));  print("")  }
#########
Counts.ORG= as.data.frame(Counts.ORG, stringsAsFactors = FALSE);
Annotation.Counts= Annotation.ORG
##################
ERROR.Targets= Targets.Full
ERROR.Targets[,"ID..Counts"]= colnames(Counts.ORG);     ERROR.Targets[,"ID..targets"]= Targets.Full[,TXT_RNAseq..COLUMN.Samples.ID]  
ERROR.Targets= ERROR.Targets[ERROR.Targets[,"ID..Counts"]!=ERROR.Targets[,"ID..targets"],];
if(nrow(ERROR.Targets)>0) {
	print("");   print(paste(" ERROR ERROR --> samples in Input.Counts not match with Input.Samples.IDs -->", nrow(ERROR.Targets) ));  
	print(ERROR.Targets[,c("ID..Counts","ID..targets")]); print("");  stop("kaOs is stopping")
}

#################################### RNAseq..Special.Genes..File.Name {
####################################
RNAseq..Data.Special.Genes.list= c();
if(length(RNAseq..Special.Genes..File.Name)>0 | length(RNAseq..Special.Genes..PI)>0) {
	############
	Gene.MAT.ORG= as.data.frame(array(NA,c(0,3)),stringsFactors=F); 
	if(length(RNAseq..Special.Genes..File.Name)>0){ setwd(DIR_RNAseq..Inputs.Project);   Gene.MAT.ORG= read.csv(RNAseq..Special.Genes..File.Name[1]);  COL.ORG.Gene.MAT.ORG= colnames(Gene.MAT.ORG) }
	############
	RNAseq..Special.Genes..PI..USE= c()
	if(length(RNAseq..Special.Genes..PI)>0){
		RNAseq..Special.Genes..PI..USE= RNAseq..Special.Genes..PI;      RNAseq..Special.Genes..PI..USE= RNAseq..Special.Genes..PI..USE[!is.na(RNAseq..Special.Genes..PI..USE)]
		RNAseq..Special.Genes..PI..USE= unlist(strsplit(RNAseq..Special.Genes..PI..USE, ",", fixed=TRUE));     RNAseq..Special.Genes..PI..USE=RNAseq..Special.Genes..PI..USE[RNAseq..Special.Genes..PI..USE!=""];
		if(length(RNAseq..Special.Genes..PI..USE) >0){ RNAseq..Special.Genes..PI..USE= gsub("'","",RNAseq..Special.Genes..PI..USE,fixed=TRUE);   RNAseq..Special.Genes..PI..USE= gsub('"','',RNAseq..Special.Genes..PI..USE,fixed=TRUE) }
		if(length(RNAseq..Special.Genes..PI..USE) >0) RNAseq..Special.Genes..PI..USE= unique(RNAseq..Special.Genes..PI..USE[order(RNAseq..Special.Genes..PI..USE)])
		if(length(RNAseq..Special.Genes..PI..USE)<=0) RNAseq..Special.Genes..PI..USE= "No genes"
	}
	MAX.number.of.Special.Genes= max(c(1000,nrow(Gene.MAT.ORG),length(RNAseq..Special.Genes..PI..USE))) 
	RNAseq..Data.Special.Genes.list= as.data.frame(array(NA,c(MAX.number.of.Special.Genes,(length(RNAseq..Special.Genes..PI..Column.Name)+length(RNAseq..Special.Genes..List.Names)+1))),stringsFactors=F); 
	colnames(RNAseq..Data.Special.Genes.list)= c("ROWS",RNAseq..Special.Genes..PI..Column.Name, RNAseq..Special.Genes..List.Names)    
	RNAseq..Data.Special.Genes.list[,"ROWS"]= rownames(RNAseq..Data.Special.Genes.list)= c(1:MAX.number.of.Special.Genes)
	############
	Conta.File.ll=0;  MAX.genes.ll=c();		ERROR.Not.finding.Columns= c()
	if(length(RNAseq..Special.Genes..Column.Name)>0) for(Col.name.tt in RNAseq..Special.Genes..Column.Name) {	
		##########
		Conta.File.ll= Conta.File.ll + 1;    COL.List.Names..USE= RNAseq..Special.Genes..List.Names[Conta.File.ll]
		#####
		if(length(COL.ORG.Gene.MAT.ORG[COL.ORG.Gene.MAT.ORG==Col.name.tt])<=0) { ERROR.Not.finding.Columns= c(ERROR.Not.finding.Columns, Col.name.tt)
		} else {
			Gene.List.ll= as.character(Gene.MAT.ORG[,Col.name.tt]);   Gene.List.ll= Gene.List.ll[!is.na(Gene.List.ll)];    Gene.List.ll= unlist(strsplit(Gene.List.ll, ",", fixed=TRUE))
			Gene.List.ll= Gene.List.ll[Gene.List.ll!=""];  if(length(Gene.List.ll)>0){ Gene.List.ll= gsub("'","",Gene.List.ll,fixed=TRUE);   Gene.List.ll= gsub('"','',Gene.List.ll,fixed=TRUE) }
			Gene.List.ll= unique(Gene.List.ll[order(Gene.List.ll)])
			##########
			if(length(Gene.List.ll)>0) RNAseq..Data.Special.Genes.list[c(1:length(Gene.List.ll)),COL.List.Names..USE]= Gene.List.ll
			MAX.genes.ll= c(MAX.genes.ll, length(Gene.List.ll))
	}	}
	if(length(RNAseq..Special.Genes..PI..USE)>0){
		RNAseq..Data.Special.Genes.list[c(1:length(RNAseq..Special.Genes..PI..USE)),RNAseq..Special.Genes..PI..Column.Name]= RNAseq..Special.Genes..PI..USE
		##########
		MAX.genes.ll= c(MAX.genes.ll, length(RNAseq..Special.Genes..PI..USE) )
	}
	##########
	if(length(MAX.genes.ll)>0) {
		RNAseq..Data.Special.Genes.list= RNAseq..Data.Special.Genes.list[c(1:max(MAX.genes.ll)),]
		######	
		if(print..RNAseq..script>-1){ print(paste("RNAseq..Special.Genes..File.Name --> Number of lists =", (ncol(RNAseq..Data.Special.Genes.list)-1),"    genes=", nrow(RNAseq..Data.Special.Genes.list)));  print("") } 
		if(print..RNAseq..script> 1){ print(head(RNAseq..Data.Special.Genes.list));  print(tail(RNAseq..Data.Special.Genes.list)) } 
		######
		if(length(ERROR.Not.finding.Columns)>0){ print("");  print(paste(" MAIN ERROR --> csv file do not contain some collumns -->",length(ERROR.Not.finding.Columns) ));  print(ERROR.Not.finding.Columns);  print("");  stop("kaOs is sotppping") }
	} else {
		stop(paste(" ERROR WARNING --> error --> no gene and length(RNAseq..Special.Genes..File.Name)>0 | length(RNAseq..Special.Genes..PI)>0 -->",length(RNAseq..Special.Genes..File.Name),"  ",length(RNAseq..Special.Genes..PI) ))
	}
}

############################################################################################################
############################################################################################################
if(Run..Lima.models..print_tests==1) print..RNAseq..script.Bili= print..RNAseq..script else print..RNAseq..script.Bili= 1;    print..RNAseq..script.Bili= 1
######### 
Targets.Model= FFF.RNAseq..Filter.Targets..MAKE.test(Project.RNA..run, test_RUN.Model, Targets.Full, print..RNAseq..script.Bili, RNAseq..TXT.JOIN.Samples_Ratio)
test_RUN.Model..Pipe= test_RUN.Model

###########################
COL.Tisssue.libType= COL.Tisssue.libType...FF;    rm(COL.Tisssue.libType...FF);
#########    
Targets.Model[,"libType"]= Targets.Model[,COL.Tisssue.libType];
#########
Counts.Model= Counts.ORG[,as.character(Targets.Model[,TXT_RNAseq..COLUMN.Samples.ID])];		COLUMNS.ORG..Counts.Model= colnames(Counts.Model)
#########
CONDITIONs.DESeq= unique(Targets.Model[,"Test"]);    if(length(CONDITIONs.DESeq)!=2) stop(paste("ERROR ERROR   -->  length(CONDITION)!=2     -->",length(CONDITIONs.DESeq) ))

############################################################################################################
############################################################################################################ design.matrix (Targets.Model)
##################
LIB.NN= length(unique(Targets.Model[,COL.Tisssue.libType]));  if(LIB.NN>=3) stop(paste("KaOs --> this script is only for max 2, sorry :(  -->",LIB.NN))
##################
AA= length(unique(Targets.Model[Targets.Model[,"Test"]==CONDITIONs.DESeq[1],COL.Tisssue.libType]))
BB= length(unique(Targets.Model[Targets.Model[,"Test"]==CONDITIONs.DESeq[2],COL.Tisssue.libType]))
if(length(unique(Targets.Model[,COL.Tisssue.libType]))>=2 & (AA>1 | BB>1)) USE.diff.libType.DESeq=1 else USE.diff.libType.DESeq=0
if(RNAseq..FORCE_USE.only.one.libType!=0 & USE.diff.libType.DESeq==0) Repetivitive_run..only.one.libType=1 else Repetivitive_run..only.one.libType=0
if(RNAseq..FORCE_USE.only.one.libType!=0) USE.diff.libType.DESeq= 0
##################
Experiment= factor(Targets.Model[,"Test"] );    	 libType= factor(Targets.Model[,COL.Tisssue.libType])
#########
if(USE.diff.libType.DESeq==0) pasillaDesign.DESeq= data.frame(row.names= COLUMNS.ORG..Counts.Model, Experiment= Targets.Model[,"Test"] )
if(USE.diff.libType.DESeq>=1) pasillaDesign.DESeq= data.frame(row.names= COLUMNS.ORG..Counts.Model, Experiment= Targets.Model[,"Test"], libType= Targets.Model[,"libType"])
#########
if(USE.diff.libType.DESeq==0){ design.Model= model.matrix(~Experiment);                		colnames(design.Model)= c("A","B")	 }
if(USE.diff.libType.DESeq>=1){ design.Model= model.matrix(~Experiment+libType);     		colnames(design.Model)= c("A","B","C") }
rownames(design.Model)=COLUMNS.ORG..Counts.Model;   if(print..RNAseq..script>1) print(design.Model)
#########
if(USE.diff.libType.DESeq==0) GROUP.Use= factor(paste0(Targets.Model[,"Test"]))
if(USE.diff.libType.DESeq>=1) GROUP.Use= factor(paste0(Targets.Model[,"Test"], ".", Targets.Model[,"libType"]))

############################################################################################################
############################################################################################################
if(Run..Lima.models..print_tests==1) {
if(Repetivitive_run..only.one.libType==1) {
	print("");  print(paste("KaOs note: --> This is a repetitive run  --> we are omitted it  -->", Repetivitive_run..only.one.libType ));  print("")
	print("#######################################################################################################");                      print("")
} else {

############################################################################################################
############################################################################################################  Boris filters
if(RNAseq..QC..filtering..Counts==1) {
	print(paste(" --> Running filters Counts ==>  Min rows=",RNAseq..QC..Min.Sum.Counts..Sample,"     Min sample/row ==>", RNAseq..QC..Min.Samples..per_test))
	#########
	RNAseq..QC..Min.Sum.Counts..Sample..USE= RNAseq..QC..Min.Sum.Counts..Sample;    RNAseq..QC..Min.Samples..per_test..USE = RNAseq..QC..Min.Samples..per_test
	#########
	DIM.old= dim(Counts.Model);   reads_clean = Counts.Model;   COL.SAMPLE.QC.ERROR=c()
	reads_wo_ctrl=reads_clean[,which(!colSums(reads_clean)< RNAseq..QC..Min.Sum.Counts..Sample..USE)];    if(DIM.old[2]!=ncol(reads_wo_ctrl)){ COL.SAMPLE.QC.ERROR= colnames(reads_wo_ctrl);  reads_wo_ctrl=reads_clean } 
	keep=apply(reads_wo_ctrl,1,function(x) length(x[x>=1])>=RNAseq..QC..Min.Samples..per_test..USE )
	Counts.Model = reads_wo_ctrl[keep,];     if(length(Annotation.Counts)>0) Annotation.Counts= Annotation.Counts[keep,]
	#########
	if(print..RNAseq..script>=0){ print(DIM.old);  print(dim(reads_wo_ctrl));  print(dim(Counts.Model));  print("") }
	#########
	if(length(COL.SAMPLE.QC.ERROR)>0) {
		QC.ERROR.TEMP= c();  for(SAM.j in COLUMNS.ORG..Counts.Model) if(length(COL.SAMPLE.QC.ERROR[COL.SAMPLE.QC.ERROR==SAM.j])<=0) QC.ERROR.TEMP= c(QC.ERROR.TEMP, SAM.j)
		print("");  print(paste(" ERROR  ERROR  -->  Samples have problems on QC -->", (DIM.old[2]-ncol(Counts.Model))));  print(QC.ERROR.TEMP);    print("kaOS:forcing to continue");   print("");  
		SAMPLES.QC.ERROR..ALL= c(SAMPLES.QC.ERROR..ALL, QC.ERROR.TEMP)
}	}

############################################################################################################
############################################################################################################
##DGE_0= DGEList(counts=Counts.Model);								if(print..RNAseq..script>0){ print(DGE_0$samples);   print(DGE_0$group) }
  DGE_0= DGEList(counts=Counts.Model, group=GROUP.Use);				if(print..RNAseq..script>0){ print(DGE_0$samples);   print(DGE_0$group) }

############################################################################################################ apply scale normalization to RNA-seq read counts
############################################################################################################
##DGE_0= calcNormFactors(DGE_0, method = c("TMM","RLE","upperquartile","none")[1])  ### boris script
  DGE_0= calcNormFactors(DGE_0)
########
if(print..RNAseq..script>1){ head(DGE$counts);  head(DGE$lib.size);  head(DGE$norm.factors);  head(DGE$samples);  head(DGE_0$group);  head(DGE$genes);  head(DGE$remove.zeros) }
########
dge_MAIN <<- DGE_0;

############################################################################################################
############################################################################################################  initializtion
CASE.Model.VEC=c();    COl.Needs.TABLEs.ALL.ALL= c("logCPM","LR","logFC","PValue","FDR", "Model","RawPro.gene")
TABLE00= as.data.frame(array(NA,c(0,length(COl.Needs.TABLEs.ALL.ALL))),stringsFactors=F); colnames(TABLE00)= COl.Needs.TABLEs.ALL.ALL
TABLE1=TABLE00;  TABLE2=TABLE00;  TABLE3=TABLE00;  TABLE4=TABLE00;  TABLE5=TABLE00;  TABLE6=TABLE00;  TABLE7=TABLE00;  TABLE8=TABLE00;  TABLE9=TABLE00;  TABLE10=TABLE00;  TABLE11=TABLE00;  TABLE12=TABLE00;

############################################################################################################
############################################################################################################  standard limma pipeline
if(length(Run.Model..LIMMA)>0) if(Run.Model..LIMMA>=1) {
	print(paste("#############################################################   Running Limma 1         ", date() ))
	#########
	CPM_lima = cpm (DGE_0, log=TRUE, prior.count=3)
	########
	fit    = lmFit(CPM_lima, design.Model)
	fit1   = eBayes(fit, trend=TRUE   )		  		## standard limma pipeline, using the trend=TRUE
	fit2   = treat (fit, lfc=log2(1.2))         	## give more weight to fold-changes in the gene ranking,
	########
	TABLE1 = topTable(fit1, coef=ncol(design.Model), sort.by="none", lfc=0, number=1000000);   if(print..RNAseq..script>1){ print(""); print("Results1:");  print(head(TABLE1));  print(dim(TABLE1))  }
	TABLE2 = topTreat(fit2, coef=ncol(design.Model), sort.by="none", lfc=0, number=1000000);   if(print..RNAseq..script>1){ print(""); print("Results2:");  print(head(TABLE2));  print(dim(TABLE2))  }  
	########
	TABLE1.Boris1= TABLE1[TABLE1[,"P.Value"]<=PRINT.PValue,];  TABLE1.Boris2= TABLE1[TABLE1[,"P.Value"]<=PRINT.PValue & TABLE1[,"adj.P.Val"]<=PRINT.FDR.Lima,];  print(paste(" RAw data=",nrow(TABLE1),"   pval<0.05=",nrow(TABLE1.Boris1),"   FDR=",nrow(TABLE1.Boris2) )) 
	TABLE2.Boris1= TABLE2[TABLE2[,"P.Value"]<=PRINT.PValue,];  TABLE2.Boris2= TABLE2[TABLE2[,"P.Value"]<=PRINT.PValue & TABLE2[,"adj.P.Val"]<=PRINT.FDR.Lima,];  print(paste(" RAw data=",nrow(TABLE2),"   pval<0.05=",nrow(TABLE2.Boris1),"   FDR=",nrow(TABLE2.Boris2) ))
	
	################################ 
	## When the library sizes are quite variable between samples, then the voom approach is theoretically more powerful than limma-trend
	########
	if(RNAseq..Models..VOOM==1) VOOM = voom(DGE_0 , design.Model, plot=FALSE)      ## (library sizes are quite variable between samples)
	if(RNAseq..Models..VOOM==2) VOOM = voom(Counts.Model, design.Model, plot=FALSE)
	if(RNAseq..Models..VOOM==3) VOOM = voom(Counts.Model, design.Model, plot=FALSE, normalize="quantile");		## If the data are very noisy
	########
	fit_Vo = lmFit(VOOM, design.Model)
	fit1_Vo= eBayes(fit_Vo, trend=TRUE   )		  		## standard limma pipeline, using the trend=TRUE
	fit2_Vo= treat (fit_Vo, lfc=log2(1.2))         		## give more weight to fold-changes in the gene ranking,
	########
	TABLE3 = topTable(fit1_Vo, coef=ncol(design.Model), sort.by="none", lfc=0, number=1000000);   if(print..RNAseq..script>1){ print(""); print("Results1:");  print(head(TABLE3));  print(dim(TABLE3))  }
	TABLE4 = topTreat(fit2_Vo, coef=ncol(design.Model), sort.by="none", lfc=0, number=1000000);   if(print..RNAseq..script>1){ print(""); print("Results2:");  print(head(TABLE4));  print(dim(TABLE4))  }  
	########
	TABLE3.Boris1= TABLE3[TABLE3[,"P.Value"]<=PRINT.PValue,];  TABLE3.Boris2= TABLE3[TABLE3[,"P.Value"]<=PRINT.PValue & TABLE3[,"adj.P.Val"]<=PRINT.FDR.Lima,];  print(paste(" RAw data=",nrow(TABLE3),"   pval<0.05=",nrow(TABLE3.Boris1),"   FDR=",nrow(TABLE3.Boris2) )) 
	TABLE4.Boris1= TABLE4[TABLE2[,"P.Value"]<=PRINT.PValue,];  TABLE4.Boris2= TABLE2[TABLE2[,"P.Value"]<=PRINT.PValue & TABLE4[,"adj.P.Val"]<=PRINT.FDR.Lima,];  print(paste(" RAw data=",nrow(TABLE4),"   pval<0.05=",nrow(TABLE4.Boris1),"   FDR=",nrow(TABLE4.Boris2) ))
	########
	################################ 
	if(nrow(TABLE1)>0){ TABLE1[,"RawPro.gene"]=rownames(TABLE1);  TABLE1[,"Model"]="Lima1"      } else TABLE1= TABLE00;     if(nrow(TABLE2)>0){ TABLE2[,"RawPro.gene"]=rownames(TABLE2);  TABLE2[,"Model"]="Lima2"      } else TABLE2= TABLE00
	if(nrow(TABLE3)>0){ TABLE3[,"RawPro.gene"]=rownames(TABLE3);  TABLE3[,"Model"]="Lima3-Voom" } else TABLE3= TABLE00;     if(nrow(TABLE4)>0){ TABLE4[,"RawPro.gene"]=rownames(TABLE4);  TABLE4[,"Model"]="Lima4-Voom" } else TABLE4= TABLE00
	CASE.Model.VEC= c(CASE.Model.VEC,"Lima1","Lima2","Lima3-Voom","Lima4-Voom")
	dge_Lima <<- CPM_lima;     fit_Lima1  <<- fit1;		fit_Lima2  <<- fit2;    dge_Lima_Voom <<- VOOM;     fit_Lima_Voom1  <<- fit1_Vo;     fit_Lima_Voom2  <<- fit2_Vo
	########
}

############################################################################################################
############################################################################################################  standard edgeR pipeline
if(length(Run.Model..edgeR)>0) if(Run.Model..edgeR>=1) for(RUN_edgeR_Standard_Boris in Run.edgeR..Standard.Boris..VEC) {
	################################################ 
	################################################ calculate dispersion
	if(RUN_edgeR_Standard_Boris==1) {
		print(""); print(paste("#############################################################   Running edgeR 1         ", date() ))
		DGE_edgeR= estimateDisp(DGE_0, design.Model, robust=TRUE);  			if(print..RNAseq..script>0) DGE_edgeR$common.dispersion
	} else {
		print(""); print(paste("#############################################################   Running edgeR 2         ", date() ))
        DGE_edgeR= estimateGLMCommonDisp (DGE_0    , design.Model, verbose=T);  # Boris	
   #####DGE_edgeR= estimateGLMCommonDisp (DGE_0    , design.Model)	
		DGE_edgeR= estimateGLMTagwiseDisp(DGE_edgeR, design.Model)
		DGE_edgeR= estimateGLMTrendedDisp(DGE_edgeR, design.Model)
	}
  	################################################ quasi-likelihood F-tests (is more popular)
  	################################################ 
  	fit_lrt= glmFit  (DGE_edgeR, design.Model);						if(print..RNAseq..script>0) head(fit_lrt$coefficients)  ## The quasi-likelihood method is highly recommended for differential expression analyses of bulk RNA-seq Data as it gives stricter error rate control by accounting for the uncertainty in dispersion estimation
	########
	TEST1.lrt= glmLRT(fit_lrt);  									if(print..RNAseq..script>1){ print(topTags(TEST1.lrt));   summary(decideTests(TEST1.lrt)) }
  ##TEST1.lrt= glmLRT(fit_lrt, coef=COEF.edgeR);  					if(print..RNAseq..script>1){ print(topTags(TEST1.lrt));   summary(decideTests(TEST1.lrt)) }
  ##TEST1.lrt= glmLRT(fit_lrt, contrast=c(0,-1))
   ########
	TEST1.lrt$table$FDR= p.adjust(TEST1.lrt$table$PValue, method="BH");
	
  	################################################ likelihood ratio tests: 
  	################################################
	fit_qlf= glmQLFit(DGE_edgeR, design.Model, robust=TRUE);   		if(print..RNAseq..script>0) head(fit_qlf$coefficients)  ## The likelihood ratio test can be useful in some special cases such as single cell RNA-seq and datasets with no replicates
	########
    TEST1.qlf= glmQLFTest(fit_qlf           );	   					if(print..RNAseq..script>1){ print(topTags(TEST1.qlf));   summary(decideTests(TEST1.qlf)) }
  ##if(USE.diff.libType.DESeq==0) COEF.edgeR= c(2);
  ##if(USE.diff.libType.DESeq>=1) COEF.edgeR= c(2,3);
  ##TEST1.qlf= glmQLFTest(fit_qlf ,coef=COEF.edgeR);	   			if(print..RNAseq..script>1){ print(topTags(TEST1.qlf));   summary(decideTests(TEST1.qlf)) }  	
 #########
  ##if(USE.diff.libType.DESeq==0){ MY.Contrasts = makeContrasts(B-A, levels=design.Model) };   ##con = makeContrasts(B.pregnant - B.lactate, levels=design.Model)
  ##if(USE.diff.libType.DESeq>=1){ MY.Contrasts = makeContrasts(B-A, levels=design.Model) }
  ##TEST1.qlf= glmQLFTest(fit_qlf , contrast=MY.Contrasts);  		if(print..RNAseq..script>0){ print(topTags(TEST1.qlf));   summary(decideTests(TEST1.qlf)) }
  	########
    TEST1.qlf$table$FDR= p.adjust(TEST1.qlf$table$PValue, method="BH")
  	################
  	TABLEXX= TEST1.lrt$table;  	TABLEXX.Boris1= TABLEXX[TABLEXX[,"PValue"]<=PRINT.PValue,];  TABLEXX.Boris2= TABLEXX[TABLEXX[,"PValue"]<=PRINT.PValue & TABLEXX[,"FDR"]<=PRINT.FDR.edgeR,];  print(paste(" RAw data=",nrow(TABLEXX),"   pval<0.05=",nrow(TABLEXX.Boris1),"   FDR=",nrow(TABLEXX.Boris2) ));   if(print..RNAseq..script>1){ print(""); print("Results1:");  print(head(TABLEXX));  print(dim(TABLEXX))  }
  	TABLEYY= TEST1.qlf$table;   TABLEYY.Boris1= TABLEYY[TABLEYY[,"PValue"]<=PRINT.PValue,];  TABLEYY.Boris2= TABLEYY[TABLEYY[,"PValue"]<=PRINT.PValue & TABLEYY[,"FDR"]<=PRINT.FDR.edgeR,];  print(paste(" RAw data=",nrow(TABLEYY),"   pval<0.05=",nrow(TABLEYY.Boris1),"   FDR=",nrow(TABLEYY.Boris2) ));   if(print..RNAseq..script>1){ print(""); print("Results2:");  print(head(TABLEYY));  print(dim(TABLEYY))  }  
  	#########
  	if(nrow(TABLEXX)>0) TABLEXX[,"RawPro.gene"]=rownames(TABLEXX) else TABLEXX= TABLE00;			if(nrow(TABLEYY)>0) TABLEYY[,"RawPro.gene"]=rownames(TABLEYY) else TABLEYY= TABLE00
	################
	if(RUN_edgeR_Standard_Boris==1) {
		#########
		if(nrow(TABLEXX)>0) TABLEXX[,"Model"]="edgeR.ltr1";     if(nrow(TABLEYY)>0) TABLEYY[,"Model"]="edgeR.qlf1";     TABLE5= TABLEXX;	 TABLE6= TABLEYY;   CASE.Model.VEC= c(CASE.Model.VEC,"edgeR.ltr1","edgeR.qlf1")
		#########
		dge_edgeR   <<- DGE_edgeR;       fit_edgeR_lrt <<- fit_lrt;        fit_edgeR_qlf   <<- fit_qlf;       res_edgeR_lrt  <<- TEST1.lrt;      res_edgeR_qlf  <<- TEST1.qlf
	} else {
		if(nrow(TABLEXX)>0) TABLEXX[,"Model"]="edgeR.ltr2";     if(nrow(TABLEYY)>0) TABLEYY[,"Model"]="edgeR.qlf2";     TABLE7= TABLEXX;	 TABLE8= TABLEYY;   CASE.Model.VEC= c(CASE.Model.VEC,"edgeR.ltr2","edgeR.qlf2")
		#########
		dge_edgeR_B <<- DGE_edgeR;       fit_edgeR_B_lrt <<- fit_lrt;      fit_edgeR_B_qlf <<- fit_qlf;       res_edgeR_B_lrt <<- TEST1.lrt;    res_edgeR_B_qlf <<- TEST1.qlf
	}
	################################################ plotting
  	################################################
  ##points <- c(0,1,2,15,16,17);  colors <- rep(c("blue", "darkgreen", "red"), 2);  
  ##plotMDS(y, col=colors[group], pch=points[group]);   legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2);  plotMDS(y, col=rep(1:2, each=3));
  ##barplot(y$samples$lib.size*1e-6, names=1:69, ylab="Library size (millions)")
  ##barcodeplot(TEST1.lrt.qlf$table$logFC, ind[[1]], main=names(ind)[1])
  ##barcodeplot(TEST1.lrt.qlf$table$logFC, index[[1]], index[[2]])
  	################
  	######## 
	#### OOO <- order(TEST1.lrt$table$PValue);  cpm(DGE_edgeR)[OOO[1:10],]
	#### OOO <- order(      lrt$table$PValue);  cpm(DGE_edgeR)[OOO[1:10],]
	#### setwd(paste(MAIN.MAIN.DIR...OUTPUT..RNseq,"XX..TEST/",sep=""));   pdf("TEST.pdf"); plotMD(TEST1.lrt); abline(h=c(-1, 1), col="blue");  plotMD(TEST1.qlf); abline(h=c(-1, 1), col="blue");  dev.off();
	########
}

############################################################################################################
############################################################################################################  Run.Models..ORG.BORIS
if(Run.Models..ORG.BORIS==1) {
	print(""); print(paste("#############################################################   Running Boris        ", date() ))
	#########
	COLUMNS.Genes.Anno= c("gene.Raw");
	COl.Format.Tables  = c("RR.table","RR.Total","TEST","logFC","logCPM","LR","FDR","AveExpr","t","P.Value","adj.P.Val","ERROR.glm","B","##Gen",COLUMNS.Genes.Anno,"gene.Synonyms","RawPro.gene","##Raw")
	#########
	design.Model2 = model.matrix(~ Test, data = Targets.Model);
	#########
	dge= DGEList(counts=Counts.Model)
	dge= calcNormFactors(dge, method = c("TMM","RLE","upperquartile","none")[1])  ### boris script
	YES= estimateGLMCommonDisp (dge, design.Model2, verbose=T)
	YES= estimateGLMTagwiseDisp(YES, design.Model2)
	#########
	fit1= glmFit(YES, design.Model2)
	test_LRT=glmLRT(fit1, contrast=c(0,-1))
	test_LRT$table$FDR = p.adjust(test_LRT$table$PValue, method="BH")
	#########
	Table.7T= test_LRT$table;     Table.7T[,"P.Value"]= Table.7T[,"PValue"];   if(nrow(Table.7T)>0){ Table.7T[,"logFC"]= -Table.7T[,"logFC"]; Table.7T[,"logCPM"]= -Table.7T[,"logCPM"]  }   ##    rr=c("84069","26155","9383")[3]; Counts.Model[rr,];  print(""); TABLE8[rr,];  Table.7T[rr,];  TABLE6[rr,];  TABLE5[rr,]   
	Table.8T= test_LRT$table;     Table.8T[,"P.Value"]= Table.8T[,"PValue"];   if(nrow(Table.8T)>0){ Table.8T[,"logFC"]= -Table.8T[,"logFC"]; Table.8T[,"logCPM"]= -Table.8T[,"logCPM"]  }   ##    ;stop();  head(Table.8T[Table.8T[,"P.Value"]<=0.001,]);   head(TABLE6[TABLE6[,"P.Value"]<=0.001,]);
	#########
	MIN.Fil.FDR..Boris=0.05;   MIN.log2.fold..Boris= c(0);  MIN.p.value..Boris=0.05
	Table.7T = FFF.RNAseq..Format.table.out(Table.7T, 7, COl.Format.Tables);  Table.7T = FFF.RNAseq..Filter.logFC.P.Value(Table.7T , MIN.log2.fold..Boris, MIN.p.value..Boris, -10);		  
	Table.8T = FFF.RNAseq..Format.table.out(Table.8T, 8, COl.Format.Tables);  Table.8T = FFF.RNAseq..Filter.logFC.P.Value(Table.8T , MIN.log2.fold..Boris, MIN.p.value..Boris, MIN.Fil.FDR..Boris);	
	#########
	print(paste(" TABLE  7 ==>", nrow(Table.7T ),"     Table  8 ==>",nrow(Table.8T ) ))
}

############################################################################################################
############################################################################################################  standard DESeq
METHOD..DESeq...estimateSizeFactors= METHOD..DESeq...estimateSizeFactors.Project
if(length(RNAseq..DESeq.METHOD..ONLY.Test.VEC[RNAseq..DESeq.METHOD..ONLY.Test.VEC==test_RUN.Model])>0) METHOD..DESeq...estimateSizeFactors= RNAseq..DESeq.METHOD..ONLY.USE.Method
##################	
if(length(Run.Model..DESeq)>0) if(Run.Model..DESeq>=1) {
	##################
	if(METHOD..DESeq...estimateSizeFactors=="NOT_RUN_DESeq"){
		     print(""); print(paste("#################       WARNING WARNING     #################   NOT RUN DESeq          method=",METHOD..DESeq...estimateSizeFactors,"              ", date() ))
	} else { print(""); print(paste("#############################################################   Running DESeq          method=",METHOD..DESeq...estimateSizeFactors,"              ", date() ))
	#########
	Counts..DESeq= Counts.Model;   for(SAM in colnames(Counts..DESeq)) Counts..DESeq[,SAM]= as.integer(Counts..DESeq[,SAM])
    #########			
	cds.DESeq1 = newCountDataSet(Counts..DESeq, pasillaDesign.DESeq$Experiment );         	if(print..RNAseq..script>1){ print("###### STEP ==> newCountDataSet");  }
	cds.DESeq1 = estimateSizeFactors( cds.DESeq1 );                             			if(print..RNAseq..script>1){ print("###### STEP ==> estimateSizeFactors");  if(print..RNAseq..script>1){ print(sizeFactors(cds.DESeq1));   print(head( counts(cds.DESeq1, normalized=TRUE ) ));  print("");  	}}
	if(METHOD..DESeq...estimateSizeFactors=="Default") cds.DESeq1 = estimateDispersions( cds.DESeq1 ) else cds.DESeq1 = estimateDispersions( cds.DESeq1, method = METHOD..DESeq...estimateSizeFactors);
	   			               				                                                if(print..RNAseq..script>1){ print("###### STEP ==> estimateDispersions");  if(print..RNAseq..script>1){      str( fitInfo(cds.DESeq1));   print(head(  fData(cds.DESeq1) ));                    print("");     }}
	res.DESeq1 = nbinomTest(cds.DESeq1, CONDITIONs.DESeq[1],CONDITIONs.DESeq[2]);  	        if(print..RNAseq..script>1){ print("###### STEP ==> nbinomTest"         );  if(print..RNAseq..script>1){ print(head(   res.DESeq1     ));  print("");	}}
	res.DESeq1[,"ERROR.glm"]=0;     
	###########################
	if(length(unique(pasillaDesign.DESeq$libType))>=2 & RNAseq..Exception..Project..libType==0) {  
		if(print..RNAseq..script>0){ print(paste("################################################### cds (conditions 1 & 2)          ", date() ));  if(print..RNAseq..script>1) print("");  }
		#####################
		cds.DESeq2= newCountDataSet(Counts..DESeq, pasillaDesign.DESeq );         			if(print..RNAseq..script>1){ print("###### STEP ==> newCountDataSet");  }
		cds.DESeq2= estimateSizeFactors( cds.DESeq2 );                             			if(print..RNAseq..script>1){ print("###### STEP ==> estimateSizeFactors");  if(print..RNAseq..script>1){ print(sizeFactors(cds.DESeq2));   print(head( counts(cds.DESeq2, normalized=TRUE ) ));  print("");  	}}
		if(METHOD..DESeq...estimateSizeFactors=="Default") cds.DESeq2= estimateDispersions( cds.DESeq2 ) else cds.DESeq2 = estimateDispersions( cds.DESeq2, method = METHOD..DESeq...estimateSizeFactors);
		                                                                                    if(print..RNAseq..script>1){ print("###### STEP ==> estimateDispersions");  if(print..RNAseq..script>1){      str( fitInfo(cds.DESeq2));   print(head(  fData(cds.DESeq2) ));  	}}
		#########
		cds.fit1.DESeq2 = fitNbinomGLMs( cds.DESeq2, count ~ libType + Experiment );       	if(print..RNAseq..script>1){ print("###### STEP ==> fitNbinomGLMs fit1");  if(print..RNAseq..script>1){ str(cds.fit1.DESeq2);    print(head(cds.fit1.DESeq2)); print("");    }}
		cds.fit0.DESeq2 = fitNbinomGLMs( cds.DESeq2, count ~ libType );						if(print..RNAseq..script>1){ print("###### STEP ==> fitNbinomGLMs fit2");  if(print..RNAseq..script>1){ str(cds.fit0.DESeq2);    print(head(cds.fit0.DESeq2)); print("");    }}
		pvalsGLM.DESeq2 = nbinomGLMTest( cds.fit1.DESeq2, cds.fit0.DESeq2 );       			if(print..RNAseq..script>1){ print("###### STEP ==> nbinomGLMTest"     );  if(print..RNAseq..script>1){ print(head(pvalsGLM.DESeq2)); print("");    }}
		padjGLM.DESeq2  = p.adjust( pvalsGLM.DESeq2, method=c("BH","bonferroni")[1]);   	if(print..RNAseq..script>1){ print("###### STEP ==> p.adjust"          );  if(print..RNAseq..script>1){ print(head(padjGLM.DESeq2 )); print("");    }} ## p.adjust to adjust for multiple testing using the Benjamini-Hochberg (BH) method.
		#####################
		res.DESeq2=res.DESeq1
		if(nrow(res.DESeq2)!=length(pvalsGLM.DESeq2) | nrow(res.DESeq2)!=length(padjGLM.DESeq2) ) stop(paste(" MAIN ERROR   -->   nrow(res.DESeq2)!=length(pvalsGLM)  --> ",nrow(res.DESeq2),"!=",length(pvalsGLM.DESeq2),"!=",length(padjGLM.DESeq2) ))
		############
		res.DESeq2[!is.na(pvalsGLM.DESeq2),"pval"]= pvalsGLM.DESeq2[!is.na(pvalsGLM.DESeq2)];    res.DESeq2[!is.na(padjGLM.DESeq2),"padj"]= padjGLM.DESeq2[!is.na(padjGLM.DESeq2)];
		############
		res.DESeq2[is.na(pvalsGLM.DESeq2),"ERROR.glm"]=1;         res.DESeq2[is.na(padjGLM.DESeq2),"ERROR.glm"]= res.DESeq2[is.na(padjGLM.DESeq2),"ERROR.glm"] + 1
		ERRROR.not.conver.DESeq2= res.DESeq2[res.DESeq2[,"ERROR.glm"]!=0,];  if(nrow(ERRROR.not.conver.DESeq2)>0) print(paste(" Warning -->   glm.fit: algorithm did not converge  --> ", nrow(ERRROR.not.conver.DESeq2) )) 
		RUN.DESeq.1b.x= 1   
	} else {
		RUN.DESeq.1b.x= 0;  cds.DESeq2= cds.DESeq1[0,];    cds.fit1.DESeq2=c();    res.DESeq2= res.DESeq1[0,];   
	}
	#########
	TABLE9 = res.DESeq1;  TABLE9.Boris1 = TABLE9 [TABLE9 [,"pval"]<=PRINT.PValue,];  TABLE9.Boris2 = TABLE9 [TABLE9 [,"pval"]<=PRINT.PValue & TABLE9 [,"padj"]<=PRINT.FDR.DESeq,];  print(paste(" RAw data=",nrow(TABLE9 ),"   pval<0.05=",nrow(TABLE9.Boris1 ),"   FDR=",nrow(TABLE9.Boris2 ) ));;   if(print..RNAseq..script>1){ print(""); print("Results1:");  print(head(TABLE9 ));  print(dim(TABLE9 ))  }
	TABLE10= res.DESeq2;  TABLE10.Boris1= TABLE10[TABLE10[,"pval"]<=PRINT.PValue,];  TABLE10.Boris2= TABLE10[TABLE10[,"pval"]<=PRINT.PValue & TABLE10[,"padj"]<=PRINT.FDR.DESeq,];  print(paste(" RAw data=",nrow(TABLE10),"   pval<0.05=",nrow(TABLE10.Boris1),"   FDR=",nrow(TABLE10.Boris2) ));;   if(print..RNAseq..script>1){ print(""); print("Results2:");  print(head(TABLE10));  print(dim(TABLE10))  }
	########
	if(nrow(TABLE9 )>0){ rownames(TABLE9 )=TABLE9 [,"id"];   TABLE9 [,"RawPro.gene"]=TABLE9 [,"id"];  TABLE9 [,"Model"]="DESeq1"  } else TABLE9 = TABLE00
	if(nrow(TABLE10)>0){ rownames(TABLE10)=TABLE10[,"id"];   TABLE10[,"RawPro.gene"]=TABLE10[,"id"];  TABLE10[,"Model"]="DESeq1b" } else TABLE10= TABLE00
	########
	CASE.Model.VEC= c(CASE.Model.VEC,"DESeq1" );  if(RUN.DESeq.1b.x==1 | 55==55) CASE.Model.VEC= c(CASE.Model.VEC,"DESeq1b")
	########
	dge_DESeq1 <<- cds.DESeq1;    dge_DESeq1b <<- cds.DESeq2;     fit_DESeq1 <<- c();     fit_DESeq1b <<- cds.fit1.DESeq2;
}}


############################################################################################################ 
############################################################################################################  DESeq2  
if(length(Run.Model..DESeq2)>0) if(Run.Model..DESeq2>=1) {
	print(paste("#############################################################   DESeq2         ", date() )) 
	#########
	DEsk2.condition..Pipe= paste("condition_",paste(rev(RNA.SEQ.Test.VEC..Pipe),collapse="_vs_"),sep="")
	#########
	Targets.DESeq2= Targets.Model[,c("Patient_number","Tissue.Column","Test")]
    Targets.DESeq2$condition <- factor(Targets.DESeq2[,"Test"])
    rownames(Targets.DESeq2) <- Targets.DESeq2[,"Patient_number"]
    #########
    Counts..DESeq= Counts.Model;   for(SAM in colnames(Counts..DESeq)) Counts..DESeq[,SAM]= as.integer(Counts..DESeq[,SAM])
    Counts..DESeq= Counts..DESeq[,Targets.DESeq2[,"Patient_number"]];     if(ncol(Counts..DESeq)!=nrow(Targets.DESeq2)) stop("MAIN ERROR ---> ncol(Counts..DESeq)!=nrow(Targets.DESeq2)")
    #########
    dds <- DESeqDataSetFromMatrix(Counts..DESeq, Targets.DESeq2, design = ~ condition)
	if(METHOD..DESeq2..mean.dispersion.Project=="Default") dds <- DESeq(dds) else dds <- DESeq(dds, fitType=METHOD..DESeq2..mean.dispersion.Project)
    CONDICTION.YES= resultsNames(dds)   ## lists the coefficients
    #########
    CONDICTION.YES= CONDICTION.YES[grep("condition_", CONDICTION.YES,fixed=TRUE)];  if(length(CONDICTION.YES)!=1) stop("MAIN ERROR --> length(CONDICTION.YES)!=1")
    DEsk2.condition..USE= DEsk2.condition..Pipe;  if( CONDICTION.YES!=DEsk2.condition..Pipe) DEsk2.condition..USE= paste("condition_",paste((RNA.SEQ.Test.VEC..Pipe),collapse="_vs_"),sep="")
    if(CONDICTION.YES!=DEsk2.condition..USE) stop("MAIN ERROR ---> CONDICTION.YES!=DEsk2.condition..USE")
    res1 <- results  (dds, name=DEsk2.condition..USE)
    res2 <- lfcShrink(dds, coef=DEsk2.condition..USE, type="apeglm")
    #########
    TABLE11= data.frame(res1);  if(nrow(TABLE11)>0){ TABLE11[,"RawPro.gene"]=rownames(TABLE11);  TABLE11[,"Model"]="DESeq2"   } else TABLE11= TABLE00 
    TABLE12= data.frame(res2);  if(nrow(TABLE12)>0){ TABLE12[,"RawPro.gene"]=rownames(TABLE12);  TABLE12[,"Model"]="DESeq2b"  } else TABLE12= TABLE00
    ########
    Pvalue.NA.x= 0
  ##TABLE11= data.frame(res1);  if(nrow(TABLE11)>0){ TABLE11[is.na(TABLE11[,"pvalue"]),"pvalue"]= Pvalue.NA.x;  TABLE11[is.na(TABLE11[,"padj"]),"padj"]= Pvalue.NA.x;  TABLE11[,"RawPro.gene"]=rownames(TABLE11);  TABLE11[,"Model"]="DESeq2"   } else TABLE11= TABLE00
  ##TABLE12= data.frame(res2);  if(nrow(TABLE12)>0){ TABLE12[is.na(TABLE12[,"pvalue"]),"pvalue"]= Pvalue.NA.x;  TABLE12[is.na(TABLE12[,"padj"]),"padj"]= Pvalue.NA.x;  TABLE12[,"RawPro.gene"]=rownames(TABLE12);  TABLE12[,"Model"]="DESeq2b"  } else TABLE12= TABLE00
    ########
    TABLE11.Boris1= TABLE11.Boris2= TABLE11;   TABLE11.Boris1[is.na(TABLE11.Boris1[,"pvalue"]),"pvalue"]= Pvalue.NA.x;  TABLE11.Boris1[is.na(TABLE11.Boris1[,"padj"]),"padj"]= Pvalue.NA.x;   TABLE11.Boris2[is.na(TABLE11.Boris2[,"pvalue"]),"pvalue"]= Pvalue.NA.x;  TABLE11.Boris2[is.na(TABLE11.Boris2[,"padj"]),"padj"]= Pvalue.NA.x
    TABLE12.Boris1= TABLE12.Boris2= TABLE12;   TABLE12.Boris1[is.na(TABLE12.Boris1[,"pvalue"]),"pvalue"]= Pvalue.NA.x;  TABLE12.Boris1[is.na(TABLE12.Boris1[,"padj"]),"padj"]= Pvalue.NA.x;   TABLE12.Boris2[is.na(TABLE12.Boris2[,"pvalue"]),"pvalue"]= Pvalue.NA.x;  TABLE12.Boris2[is.na(TABLE12.Boris2[,"padj"]),"padj"]= Pvalue.NA.x
    ########
    TABLE11.Boris1= TABLE11.Boris1[TABLE11.Boris1[,"pvalue"]<=PRINT.PValue,];  TABLE11.Boris2= TABLE11.Boris2[TABLE11.Boris2[,"pvalue"]<=PRINT.PValue & TABLE11.Boris2[,"padj"]<=PRINT.FDR.DESeq,];  print(paste(" RAw data=",nrow(TABLE11),"   pval<0.05=",nrow(TABLE11.Boris1),"   FDR=",nrow(TABLE11.Boris2) ));;   if(print..RNAseq..script>1){ print(""); print("Results1:");  print(head(TABLE11));  print(dim(TABLE11))  }
    TABLE12.Boris1= TABLE12.Boris1[TABLE12.Boris1[,"pvalue"]<=PRINT.PValue,];  TABLE12.Boris2= TABLE12.Boris2[TABLE12.Boris2[,"pvalue"]<=PRINT.PValue & TABLE12.Boris2[,"padj"]<=PRINT.FDR.DESeq,];  print(paste(" RAw data=",nrow(TABLE12),"   pval<0.05=",nrow(TABLE12.Boris1),"   FDR=",nrow(TABLE12.Boris2) ));;   if(print..RNAseq..script>1){ print(""); print("Results2:");  print(head(TABLE12));  print(dim(TABLE12))  }
	#################
	if(METHOD..DESeq2..mean.dispersion.Project=="Default") Desq2..vsd= vst(dds, blind=FALSE) else Desq2..vsd= vst(dds, blind=FALSE, fitType=METHOD..DESeq2..mean.dispersion.Project) 
    dge_DESeq2 <<- dds;  Desq2..vsd <<- Desq2..vsd;  CASE.Model.VEC= c(CASE.Model.VEC,"DESeq2" ,"DESeq2b") 
}

############################################################################################################ CASE.Model.VEC 
############################################################################################################
print(""); print(paste("#############################################################   CASE.Model.VEC       ", date() ))
##########
if(length(RNAseq..Models..NOT_RUN[RNAseq..Models..NOT_RUN=="Lima" ])>0) for(RM.mod.x in RNAseq..Models.Alg..Lima.VEC ) RNAseq..Models.Alg..Save.Pipe= RNAseq..Models.Alg..Save.Pipe[RNAseq..Models.Alg..Save.Pipe!=RM.mod.x]
if(length(RNAseq..Models..NOT_RUN[RNAseq..Models..NOT_RUN=="edgeR"])>0) for(RM.mod.x in RNAseq..Models.Alg..edgeR.VEC) RNAseq..Models.Alg..Save.Pipe= RNAseq..Models.Alg..Save.Pipe[RNAseq..Models.Alg..Save.Pipe!=RM.mod.x]
if(length(RNAseq..Models..NOT_RUN[RNAseq..Models..NOT_RUN=="DESeq"])>0) for(RM.mod.x in RNAseq..Models.Alg..DESeq.VEC) RNAseq..Models.Alg..Save.Pipe= RNAseq..Models.Alg..Save.Pipe[RNAseq..Models.Alg..Save.Pipe!=RM.mod.x]
##########
if(length(RNAseq..Models.Alg..Save.Pipe)<=0) stop("MAIN ERROR ----> length(RNAseq..Models.Alg..Save.Pipe)<=0") else RNAseq..Models.Alg..Save.Pipe= unique(RNAseq..Models.Alg..Save.Pipe)
##########
ERRRR.notfind=c();   		for(MOD.xx in RNAseq..Models.Alg..Save.Pipe) if(length(CASE.Model.VEC               [CASE.Model.VEC               ==MOD.xx])<=0) ERRRR.notfind     =  c(ERRRR.notfind      , MOD.xx)
Models.remove..Pipe=c();   	for(MOD.xx in CASE.Model.VEC               ) if(length(RNAseq..Models.Alg..Save.Pipe[RNAseq..Models.Alg..Save.Pipe==MOD.xx])<=0) Models.remove..Pipe= c(Models.remove..Pipe, MOD.xx)
####
if(length(RNAseq..Models..NOT_RUN[RNAseq..Models..NOT_RUN=="Lima" ])>0) Models.remove..Pipe=c(Models.remove..Pipe, RNAseq..Models.Alg..Lima.VEC )
if(length(RNAseq..Models..NOT_RUN[RNAseq..Models..NOT_RUN=="edgeR"])>0) Models.remove..Pipe=c(Models.remove..Pipe, RNAseq..Models.Alg..edgeR.VEC)
if(length(RNAseq..Models..NOT_RUN[RNAseq..Models..NOT_RUN=="DESeq"])>0) Models.remove..Pipe=c(Models.remove..Pipe, RNAseq..Models.Alg..DESeq.VEC);   if(length(Models.remove..Pipe)>0) Models.remove..Pipe=unique(Models.remove..Pipe)
##########
if(length(ERRRR.notfind      )>0){  print("");  print(paste(" MAIN ERROR ----> models not find ---> length(ERRRR.notfind)>0  -->"  , length(ERRRR.notfind) ));  print(ERRRR.notfind);  print("");  print("");  stop("kaOS is stopping") }
if(length(Models.remove..Pipe)>0){  print("");  print(paste(" WARNING WARNING ----> models not consider in saving Max,Mean,Min -->", length(Models.remove..Pipe)          ,"     -->", paste(Models.remove..Pipe          , collapse=",") ));   }
                                    print("");  print(paste("                 ----> models yes consider in saving Max,Mean,Min -->", length(RNAseq..Models.Alg..Save.Pipe),"     -->", paste(RNAseq..Models.Alg..Save.Pipe, collapse=",") ));

############################################################################################################ TABLEs.ALL.ALL 
############################################################################################################
if(nrow(TABLE1 )>0){ TABLE1 [,"FDR"]= TABLE1 [,"adj.P.Val"];     TABLE1 [,"PValue"]= TABLE1 [,"P.Value"];     TABLE1 [,"logCPM"]= -10; 		TABLE1 [,"LR" ]= -10  }      
if(nrow(TABLE2 )>0){ TABLE2 [,"FDR"]= TABLE2 [,"adj.P.Val"];     TABLE2 [,"PValue"]= TABLE2 [,"P.Value"];     TABLE2 [,"logCPM"]= -10; 		TABLE2 [,"LR" ]= -10  }
if(nrow(TABLE3 )>0){ TABLE3 [,"FDR"]= TABLE3 [,"adj.P.Val"];     TABLE3 [,"PValue"]= TABLE3 [,"P.Value"];     TABLE3 [,"logCPM"]= -10; 		TABLE3 [,"LR" ]= -10  }      
if(nrow(TABLE4 )>0){ TABLE4 [,"FDR"]= TABLE4 [,"adj.P.Val"];     TABLE4 [,"PValue"]= TABLE4 [,"P.Value"];     TABLE4 [,"logCPM"]= -10; 		TABLE4 [,"LR" ]= -10  }
if(nrow(TABLE5 )>0){  a=1                                                                                                                                         } ## TABLE5[,"logFC"]= -TABLE5[,"logFC"]; TABLE5[,"logCPM"]= -TABLE5[,"logCPM"]
if(nrow(TABLE6 )>0){                                                                                                                   		TABLE6 [,"LR" ]= -10  } ## TABLE7[,"logFC"]= -TABLE7[,"logFC"]; TABLE7[,"logCPM"]= -TABLE7[,"logCPM"]
if(nrow(TABLE7 )>0){  a=1 }
if(nrow(TABLE8 )>0){                                                                                                                   		TABLE8 [,"LR" ]= -10  } 
if(nrow(TABLE9 )>0){ TABLE9 [,"FDR"]= TABLE9 [,"padj"];          TABLE9 [,"PValue"]= TABLE9 [,"pval"];        TABLE9 [,"logCPM"]= -10;  	TABLE9 [,"LR" ]= -10;  	TABLE9 [,"logFC"]= TABLE9 [,"log2FoldChange"]    }
if(nrow(TABLE10)>0){ TABLE10[,"FDR"]= TABLE10[,"padj"];          TABLE10[,"PValue"]= TABLE10[,"pval"];        TABLE10[,"logCPM"]= -10;  	TABLE10[,"LR" ]= -10; 	TABLE10[,"logFC"]= TABLE10[,"log2FoldChange"]    } 
if(nrow(TABLE11)>0){ TABLE11[,"FDR"]= TABLE11[,"padj"];          TABLE11[,"PValue"]= TABLE11[,"pvalue"];      TABLE11[,"logCPM"]= -10;  	TABLE11[,"LR" ]= -10;  	TABLE11[,"logFC"]= TABLE11[,"log2FoldChange"]    }
if(nrow(TABLE12)>0){ TABLE12[,"FDR"]= TABLE12[,"padj"];          TABLE12[,"PValue"]= TABLE12[,"pvalue"];      TABLE12[,"logCPM"]= -10;  	TABLE12[,"LR" ]= -10; 	TABLE12[,"logFC"]= TABLE12[,"log2FoldChange"]    }  
############
COl.Needs= COl.Needs.TABLEs.ALL.ALL
TABLEs.ALL.ALL  = rbind(    TABLE1[,COl.Needs],TABLE2[,COl.Needs],TABLE3[,COl.Needs],TABLE4[,COl.Needs],TABLE5[,COl.Needs],TABLE6[,COl.Needs],TABLE7[,COl.Needs],TABLE8[,COl.Needs],TABLE9[,COl.Needs],TABLE10[,COl.Needs],TABLE11[,COl.Needs],TABLE12[,COl.Needs])
MAX.NUM_of_ROWS= max(c(nrow(TABLE1),      nrow(TABLE2),      nrow(TABLE3),      nrow(TABLE4),      nrow(TABLE5),      nrow(TABLE6),      nrow(TABLE7),      nrow(TABLE8),      nrow(TABLE9),      nrow(TABLE10),      nrow(TABLE11),      nrow(TABLE12)          ));
############
if(nrow(TABLEs.ALL.ALL)>0) TABLEs.ALL.ALL[,"Test.XX"]= test_RUN.Model;
############
############ to make AFF supress-induced
if(RNAseq..USE.Reverse.Targets.Model==0) if(nrow(TABLEs.ALL.ALL)>0) TABLEs.ALL.ALL[,"logFC"]= -TABLEs.ALL.ALL[,"logFC"]

############################################################################################################ Join tables (TABLE.ALL.JOIN) {
############################################################################################################
print("");    print(paste("#############################################################   Join models         ", date() ))

################################################ 
print..Res.JOIN= c(1,0)[1];
RUN.Togo= c(1,2,11,12)[1];  TABLE.ALL.JOIN = FFF.RNAseq..JOIN.Algorithms.Tables.Togo(RUN.Togo, TABLEs.ALL.ALL, test_RUN.Model, ".Max", -1, -1, -1, print..Res.JOIN, c(".Mean",".Med",".Min"))
RUN.Togo= c(1,2,11,12)[3];  TABLE.ALL.JOIN1= FFF.RNAseq..JOIN.Algorithms.Tables.Togo(RUN.Togo, TABLEs.ALL.ALL, test_RUN.Model, ".Max", -1, -1, -1, print..Res.JOIN)
RUN.Togo= c(1,2,11,12)[4];  TABLE.ALL.JOIN2= FFF.RNAseq..JOIN.Algorithms.Tables.Togo(RUN.Togo, TABLEs.ALL.ALL, test_RUN.Model, ".Max", -1, -1, -1, print..Res.JOIN)
rm(TXT.Main.filters..FF, COL.SAVE..TABLE..pvlues..FF, COL.SAVE..TABLE.test.NEW..FF) 
############
if(MAX.NUM_of_ROWS!=nrow(TABLE.ALL.JOIN) | MAX.NUM_of_ROWS!=nrow(TABLE.ALL.JOIN1) | MAX.NUM_of_ROWS!=nrow(TABLE.ALL.JOIN2) ) stop(paste("ERROR -- ERROR --> MAX.NUM_of_ROWS!=TABLE.ALL.JOIN -->", MAX.NUM_of_ROWS,"!=",nrow(TABLE.ALL.JOIN),"!=",nrow(TABLE.ALL.JOIN1),"!=",nrow(TABLE.ALL.JOIN2) ))

################################################ Table..Annotation
################################################ 
print.Anno=c(1,0)[1];       Run..all.extras= c(1,2)[1];		  Extra.columns.SAVE= c("gene.Raw","gene.Ens","gene.Synonyms","Raw.geneID");	 if(RNAseq..XXX.Gene.ID.Version=="ensembl") Extra.columns.SAVE= c(Extra.columns.SAVE,"gene_biotype")
TABLE.ALL.JOIN= FFF.RNAseq..Table..Annotation(TABLE.ALL.JOIN, Annotation.ORG, "RawPro.gene", "RawPro.gene", Run..all.extras, Extra.columns.SAVE, print.Anno)
for(COl.G in c("gene.Raw","gene.Ens","gene.Synonyms")[2:3]) TABLE.ALL.JOIN[is.na(TABLE.ALL.JOIN[,COl.G]),COl.G]=""
XXX.COL.TEMP= colnames(TABLE.ALL.JOIN);  XXX.COL.TEMP[XXX.COL.TEMP==paste("##T",test_RUN.Model,sep="")]= "##logFC";  colnames(TABLE.ALL.JOIN)= XXX.COL.TEMP 
COL.SAVE..TABLE.ALL= colnames(TABLE.ALL.JOIN);

################################################ Pvalue, ##FDR
################################################
COL.Algoritmns= c(paste("##T",test_RUN.Model,sep=""),paste("T",test_RUN.Model,".",c(RNAseq..Models.Alg..Lima.VEC, RNAseq..Models.Alg..edgeR.VEC, RNAseq..Models.Alg..DESeq.VEC),sep=""))
TABLE.ALL.JOIN1= TABLE.ALL.JOIN1[rownames(TABLE.ALL.JOIN),COL.Algoritmns];   colnames(TABLE.ALL.JOIN1)= paste(COL.Algoritmns," " ,sep="");  colnames(TABLE.ALL.JOIN1)[1]= "##Pvalue"
TABLE.ALL.JOIN2= TABLE.ALL.JOIN2[rownames(TABLE.ALL.JOIN),COL.Algoritmns];   colnames(TABLE.ALL.JOIN2)= paste(COL.Algoritmns,"  ",sep="");  colnames(TABLE.ALL.JOIN2)[1]= "##FDR"
#######
TABLE.ALL.JOIN= cbind(TABLE.ALL.JOIN, TABLE.ALL.JOIN1, TABLE.ALL.JOIN2);  
COL.SAVE..TABLE.ALL= c(COL.SAVE..TABLE.ALL[COL.SAVE..TABLE.ALL!="##End"], colnames(TABLE.ALL.JOIN1), colnames(TABLE.ALL.JOIN2),"##End") 
#######
ROWS.ORG= TABLE.ALL.JOIN[,"RawPro.gene"];  ERRO.JOIN.KK= c(ROWS.ORG[ROWS.ORG!=rownames(TABLE.ALL.JOIN1)], ROWS.ORG[ROWS.ORG!=rownames(TABLE.ALL.JOIN2)])
if(length(ERRO.JOIN.KK)>0) stop("ERROR --- > match JOIN 1 or 2 does not match");   rm(TABLE.ALL.JOIN1,TABLE.ALL.JOIN2)

################################################ Running  Mean columns {
################################################
COL.Algoritmns..FORMAT  = paste("T",test_RUN.Model,".",RNAseq..Models.Alg..Save.Pipe   ,sep="");
COL.Algoritmns..MAX.Mean= paste("T",test_RUN.Model,".",RNAseq..Models.Alg..Max.Min.Mean,sep="");    print("");
#### 
if(nrow(TABLE.ALL.JOIN)>0) for(CC.mean in c("logFC","PValue","FDR")){
	####
	if(CC.mean=="logFC" ){ COL.Algoritmns.USE= paste(COL.Algoritmns..MAX.Mean,""  ,sep="");    COL.Algoritmns.USE.FORMAT= paste(COL.Algoritmns..FORMAT,""  ,sep="");   } 	
	if(CC.mean=="PValue"){ COL.Algoritmns.USE= paste(COL.Algoritmns..MAX.Mean," " ,sep="");    COL.Algoritmns.USE.FORMAT= paste(COL.Algoritmns..FORMAT," " ,sep="");   }
	if(CC.mean=="FDR"   ){ COL.Algoritmns.USE= paste(COL.Algoritmns..MAX.Mean,"  ",sep="");    COL.Algoritmns.USE.FORMAT= paste(COL.Algoritmns..FORMAT,"  ",sep="");   }   	
	####
	print(paste(" ---> Running  Mean columns (",CC.mean,")     rows=",nrow(TABLE.ALL.JOIN),  date() ))
	####
	for(RR.mean in c(1:nrow(TABLE.ALL.JOIN))){
		MEAN.OUT.JJ=NA;    MAX.OUT.JJ=NA;    MIN.OUT.JJ=NA;    MEDIAN.OUT.JJ=NA;    TEMP.JJ= unname(unlist(c(TABLE.ALL.JOIN[RR.mean,COL.Algoritmns.USE])));    TEMP.JJ= TEMP.JJ[!is.na(TEMP.JJ)];  TEMP.JJ= TEMP.JJ[TEMP.JJ!="." & TEMP.JJ!=""];
		if(length(TEMP.JJ)>0){ TEMP.JJ= as.numeric(TEMP.JJ);      TEMP.JJ= TEMP.JJ[TEMP.JJ!=Inf];   TEMP.JJ= TEMP.JJ[TEMP.JJ!=-Inf] }
		if(length(TEMP.JJ)>0){ MEAN.OUT.JJ= mean(TEMP.JJ);        MEDIAN.OUT.JJ= median(TEMP.JJ);    if(CC.mean=="PValue" |  CC.mean=="FDR" | MEDIAN.OUT.JJ<0 )  MIN.OUT.JJ = max(TEMP.JJ) else MIN.OUT.JJ = min(TEMP.JJ) 
								                                                                     if(CC.mean=="PValue" |  CC.mean=="FDR" | MEDIAN.OUT.JJ<0 )  MAX.OUT.JJ = min(TEMP.JJ) else MAX.OUT.JJ = max(TEMP.JJ) }
		TABLE.ALL.JOIN[RR.mean,paste(CC.mean,".Mean",sep="")]= MEAN.OUT.JJ;     TABLE.ALL.JOIN[RR.mean,paste(CC.mean,".Med" ,sep="")]= MEDIAN.OUT.JJ;   
		TABLE.ALL.JOIN[RR.mean,paste(CC.mean,".Max" ,sep="")]= MAX.OUT.JJ;		TABLE.ALL.JOIN[RR.mean,paste(CC.mean,".Min" ,sep="")]= MIN.OUT.JJ	
	}
	####
	if(CC.mean=="logFC") for(CC.for in COL.Algoritmns.USE.FORMAT) TABLE.ALL.JOIN[!is.na(TABLE.ALL.JOIN[,CC.for]) & abs(TABLE.ALL.JOIN[,CC.for])!=Inf,CC.for]= as.integer(TABLE.ALL.JOIN[!is.na(TABLE.ALL.JOIN[,CC.for]) & abs(TABLE.ALL.JOIN[,CC.for])!=Inf,CC.for]*10 )/10	
	if(CC.mean!="logFC") for(CC.for in COL.Algoritmns.USE.FORMAT) TABLE.ALL.JOIN[,CC.for]= as.numeric(TABLE.ALL.JOIN[,CC.for])
	if(CC.mean=="logFC") for(CC.for in paste(CC.mean,c(".Max",".Mean",".Med",".Min"),sep="")) TABLE.ALL.JOIN[!is.na(TABLE.ALL.JOIN[,CC.for]) & abs(TABLE.ALL.JOIN[,CC.for])!=Inf,CC.for]= as.integer(TABLE.ALL.JOIN[!is.na(TABLE.ALL.JOIN[,CC.for]) & abs(TABLE.ALL.JOIN[,CC.for])!=Inf,CC.for]*100)/100
	if(CC.mean!="logFC") for(CC.for in paste(CC.mean,c(".Max",".Mean",".Med",".Min"),sep="")) TABLE.ALL.JOIN[,CC.for]= as.numeric(TABLE.ALL.JOIN[,CC.for])
}
#### 
if(nrow(TABLE.ALL.JOIN)>0){
	############
	TABLE.ALL.JOIN[,"logFC.Diff"]="";  
	TABLE.ALL.JOIN[abs(as.numeric(TABLE.ALL.JOIN[,"logFC.Max"])-as.numeric(TABLE.ALL.JOIN[,"logFC.Min" ]))>2.5,"logFC.Diff"]="Warning (Min)";
	TABLE.ALL.JOIN[abs(as.numeric(TABLE.ALL.JOIN[,"logFC.Max"])-as.numeric(TABLE.ALL.JOIN[,"logFC.Mean"]))>1.5,"logFC.Diff"]="Warning (Mean)"
	TABLE.ALL.JOIN[abs(as.numeric(TABLE.ALL.JOIN[,"logFC.Max"])-as.numeric(TABLE.ALL.JOIN[,"logFC.Med" ]))>1.5,"logFC.Diff"]="Warning (Med)"
	if(print..RNAseq..script>0){ print("");  print("Results for logFC.Diff -->");   print(table(TABLE.ALL.JOIN[,"logFC.Diff"])) }
	######
	COL.SAVE..TABLE.ALL= c(COL.SAVE..TABLE.ALL, "logFC.Diff")
}
############
if(nrow(TABLE.ALL.JOIN)>0) for(Max.Col in c(".Max",".Mean",".Med",".Min")) {
	#########
	Fileter.Col= paste( "Filters",Max.Col,sep="");      logFC.Col= paste( "logFC",Max.Col,sep="");   PValue.Col= paste( "PValue",Max.Col,sep="");      FDR.Col= paste( "FDR",Max.Col,sep="");      
	#########
	TABLE.ALL.JOIN[,Fileter.Col]= ""
	if(MAX.FDR..Filters1  >0 & MAX.FDR..Filters1  <1) TABLE.ALL.JOIN[    as.numeric(TABLE.ALL.JOIN[, FDR.Col   ])  > MAX.FDR..Filters1 & as.numeric(TABLE.ALL.JOIN[,FDR.Col])<=MAX.FDR..Filters2,Fileter.Col]= "LOW"
	if(MAX.FDR..Filters2  >0 & MAX.FDR..Filters2  <1) TABLE.ALL.JOIN[    as.numeric(TABLE.ALL.JOIN[, FDR.Col   ])  > MAX.FDR..Filters2  , Fileter.Col ]= "NOT PASS"
	if(MAX.PValue..Filters>0 & MAX.PValue..Filters<1) TABLE.ALL.JOIN[    as.numeric(TABLE.ALL.JOIN[, PValue.Col])  > MAX.PValue..Filters, Fileter.Col ]= "NOT PASS"
	if(MIN.logFC..Filters >0                        ) TABLE.ALL.JOIN[abs(as.numeric(TABLE.ALL.JOIN[, logFC.Col ])) < MIN.logFC..Filters , Fileter.Col ]= "NOT PASS"
}

################################################
################################################ RNAseq..Data.Special.Genes.list
COl.Special.gene.List= c()
if(nrow(TABLE.ALL.JOIN)>0 & length(RNAseq..Data.Special.Genes.list)>0) {
	COl.Special.gene.List= colnames(RNAseq..Data.Special.Genes.list); COl.Special.gene.List= COl.Special.gene.List[COl.Special.gene.List!="ROWS"];     TABLE.ALL.JOIN[,COl.Special.gene.List]=0;  print("");  
	##########
	for(GE.CO in COl.Special.gene.List) {
		TABLE.ALL.JOIN[,GE.CO]="";     GENE.LIST.check= RNAseq..Data.Special.Genes.list[,GE.CO];  GENE.LIST.check= GENE.LIST.check[!is.na(GENE.LIST.check)];  GENE.LIST.check= GENE.LIST.check[GENE.LIST.check!=""]
		if(length(GENE.LIST.check)>0) for(Gene.X in unique(GENE.LIST.check)) for(COl.G in c("gene.Raw","gene.Ens","gene.Synonyms")) TABLE.ALL.JOIN[!is.na(TABLE.ALL.JOIN[,COl.G]) & tolower(TABLE.ALL.JOIN[,COl.G])==tolower(Gene.X),GE.CO]="XX"
		######
		print(paste(" Special gene List --> ", GE.CO,"     genes=", length(GENE.LIST.check),"    in data=", length(TABLE.ALL.JOIN[TABLE.ALL.JOIN[,GE.CO]!="",GE.CO]) ))
}	}
#### TABLEs.ALL.ALL[TABLEs.ALL.ALL[,"gene.PLOT"]=="Sema3c",];   TABLE.ALL.JOIN[TABLE.ALL.JOIN[,"gene.PLOT"]=="Sema3c",]

################################################
################################################ COL.SAVE..TABLE.ALL
######### delete extra columns
COl..DEL.EXTRAS= c(paste("Test",test_RUN.Model,".",c("S","A","PValue","FDR","logFC"),sep=""),paste("TEST",test_RUN.Model,".",RNAseq..Models.Alg..Unique.VEC,sep=""))
for(DEL.here in c(COl..DEL.EXTRAS,"pvalue.A")) COL.SAVE..TABLE.ALL= COL.SAVE..TABLE.ALL[COL.SAVE..TABLE.ALL!=DEL.here]
#########
######### Counts.Aff","Counts.UnA
if(length(RNA.SEQ.Test.VEC..Pipe)!=2) stop("MAIN ERROR  --> length(RNA.SEQ.Test.VEC..Pipe)!=2")
COL.Counts.Aff= "Counts.Aff";  											COL.Counts.UnA= "Counts.UnA"
COL.Counts.Aff= paste("Counts.",RNA.SEQ.Test.VEC..Pipe[1],sep="");  	COL.Counts.UnA= paste("Counts.",RNA.SEQ.Test.VEC..Pipe[2],sep="")
#########
NN1.K= grep("##A",COL.SAVE..TABLE.ALL,fixed=TRUE); if(length(NN1.K)!=1) stop("ERROR grep --> ##A")
COL.SAVE..TABLE.ALL= c(COL.SAVE..TABLE.ALL[1:(NN1.K-1)],"Gene.type","##Mean",COL.Counts.Aff,COL.Counts.UnA,COL.SAVE..TABLE.ALL[(NN1.K):length(COL.SAVE..TABLE.ALL)])
#########
NN.logFC.M= grep("logFC.Max" ,COL.SAVE..TABLE.ALL,fixed=TRUE); if(length(NN.logFC.M)!=1) stop("ERROR grep --> logFC.Max" )
NN.logFC.D= grep("logFC.Diff",COL.SAVE..TABLE.ALL,fixed=TRUE); if(length(NN.logFC.D)!=1) stop("ERROR grep --> logFC.Diff")
COL.SAVE..TABLE.ALL= COL.SAVE..TABLE.ALL[unique(c( c(1:(NN.logFC.M-1)), NN.logFC.D, c((NN.logFC.M):length(COL.SAVE..TABLE.ALL) )))]
#########
if(nrow(TABLE.ALL.JOIN)>0){ TABLE.ALL.JOIN[,c("##Mean",COL.Counts.Aff,COL.Counts.UnA)]= NA;   TABLE.ALL.JOIN[,"Model"]= "Join.Alg";    TABLE.ALL.JOIN[,"Test.XX"]= test_RUN.Model }
#########
######### RNAseq..Special.Genes..PI..Column.Name, RNAseq..Special.Genes..List.Names
COL.SAVE..TABLE.ALL= c(RNAseq..Special.Genes..PI..Column.Name, RNAseq..Special.Genes..List.Names, COL.SAVE..TABLE.ALL)
#########
## head(TABLE.ALL.JOIN[,COL.SAVE..TABLE.ALL])

############################################################################################################ Annotations {
############################################################################################################
print(""); print(paste("#############################################################   Annotations         ", date() ))
############
Counts.gene.PLOT= Counts.Model;    Counts.gene.PLOT[,"RawPro.gene"]= rownames(Counts.gene.PLOT)
############
print.Anno=c(1,0)[1];  
############ 
COL.entrezgene= c();	 if(RNAseq..XXX.Gene.ID.Version=="ensembl" & DIR_Annotations.Organism!="Homo_sapiens") COL.entrezgene= "entrezgene";   
Run..all.extras= c(1,2)[2];       Extra.columns.Plotting= c("gene.Raw","gene.Ens","type_of_gene","gene_biotype",COL.entrezgene,"Raw.geneID");
TABLEs.ALL.ALL  = FFF.RNAseq..Table..Annotation(TABLEs.ALL.ALL  , Annotation.ORG, "RawPro.gene", "RawPro.gene", Run..all.extras, Extra.columns.Plotting, print.Anno)
Counts.gene.PLOT= FFF.RNAseq..Table..Annotation(Counts.gene.PLOT, Annotation.ORG, "RawPro.gene", "RawPro.gene", Run..all.extras, Extra.columns.Plotting, print.Anno)
############
TABLEs.ALL.ALL  [,"gene.PLOT"   ]= TABLEs.ALL.ALL  [,"gene.Raw"  ];       	if(RNAseq..XXX.Gene.ID.Version=="ensembl") TABLEs.ALL.ALL  [,"gene.PLOT"   ]= TABLEs.ALL.ALL  [,"gene.Ens"] 
TABLE.ALL.JOIN  [,"gene.PLOT"   ]= TABLE.ALL.JOIN  [,"gene.Raw"  ];      	if(RNAseq..XXX.Gene.ID.Version=="ensembl") TABLE.ALL.JOIN  [,"gene.PLOT"   ]= TABLE.ALL.JOIN  [,"gene.Ens"]
Counts.gene.PLOT[,"gene.PLOT"   ]= Counts.gene.PLOT[,"gene.Raw"  ];     	if(RNAseq..XXX.Gene.ID.Version=="ensembl") Counts.gene.PLOT[,"gene.PLOT"   ]= Counts.gene.PLOT[,"gene.Ens"]  
############
if(RNAseq..XXX.Gene.ID.Version=="ensembl") TABLEs.ALL.ALL[,"type.gene.PLOT"]= TABLEs.ALL.ALL[,"gene_biotype"] else TABLEs.ALL.ALL[,"type.gene.PLOT"]= TABLEs.ALL.ALL[,"type_of_gene"]; 
if(RNAseq..XXX.Gene.ID.Version=="ensembl") TABLE.ALL.JOIN[,"type.gene.PLOT"]= TABLE.ALL.JOIN[,"gene_biotype"] else TABLE.ALL.JOIN[,"type.gene.PLOT"]= TABLE.ALL.JOIN[,"type_of_gene"]; 		
############
TABLEs.ALL.ALL  [is.na(TABLEs.ALL.ALL  [,"gene.PLOT"]),"gene.PLOT"]= paste("id=",TABLEs.ALL.ALL  [is.na(TABLEs.ALL.ALL  [,"gene.PLOT"]),"RawPro.gene"])
TABLE.ALL.JOIN  [is.na(TABLE.ALL.JOIN  [,"gene.PLOT"]),"gene.PLOT"]= paste("id=",TABLE.ALL.JOIN  [is.na(TABLE.ALL.JOIN  [,"gene.PLOT"]),"RawPro.gene"])
Counts.gene.PLOT[is.na(Counts.gene.PLOT[,"gene.PLOT"]),"gene.PLOT"]= paste("id=",Counts.gene.PLOT[is.na(Counts.gene.PLOT[,"gene.PLOT"]),"RawPro.gene"])
############
TABLEs.ALL.ALL[is.na(TABLEs.ALL.ALL[,"type.gene.PLOT"]),"type.gene.PLOT"]= "";			TABLEs.ALL.ALL[,"type.gene.PLOT"]= sub("protein-coding","protein_coding",TABLEs.ALL.ALL[,"type.gene.PLOT"],fixed=TRUE)
TABLE.ALL.JOIN[is.na(TABLE.ALL.JOIN[,"type.gene.PLOT"]),"type.gene.PLOT"]= "";			TABLE.ALL.JOIN[,"type.gene.PLOT"]= sub("protein-coding","protein_coding",TABLE.ALL.JOIN[,"type.gene.PLOT"],fixed=TRUE)
TABLE.ALL.JOIN[,"Gene.type"]= TABLE.ALL.JOIN[,"type.gene.PLOT"]
#######
STOP.ERROR.NOT.unique= 0
UNI.Ge1= length(unique(TABLE.ALL.JOIN  [,"gene.PLOT"]));   if(UNI.Ge1!=nrow(TABLE.ALL.JOIN  )){ STOP.ERROR.NOT.unique=1;  print(paste("ERROR -- ERROR --> GENES names is not unique  (TABLE.ALL.JOIN  ) -->", UNI.Ge1,"!=",nrow(TABLE.ALL.JOIN  ))) }
UNI.Ge2= length(unique(Counts.gene.PLOT[,"gene.PLOT"]));   if(UNI.Ge2!=nrow(Counts.gene.PLOT)){ STOP.ERROR.NOT.unique=1;  print(paste("ERROR -- ERROR --> GENES names is not unique  (Counts.gene.PLOT) -->", UNI.Ge2,"!=",nrow(Counts.gene.PLOT))) }
if(STOP.ERROR.NOT.unique==1){ if(RNAseq..XXX.Gene.ID.Version!="ensembl") stop("kaOs is stopping") else  print("forcing to continue")  }
##########################
Counts.gene.PLOT = Counts.gene.PLOT[,c(COLUMNS.ORG..Counts.Model,"RawPro.gene","gene.PLOT")]
##########################
TABLE.ALL.JOIN.plot  = TABLE.ALL.JOIN[,c("PValue.Max","FDR.Max","logFC.Max","Filters.Max", "PValue.Mean","FDR.Mean","logFC.Mean","Filters.Mean", "PValue.Med","FDR.Med","logFC.Med","Filters.Med", "PValue.Min","FDR.Min","logFC.Min","Filters.Min", 
                                         "Model","RawPro.gene","Raw.geneID","Test.XX","gene.PLOT","type.gene.PLOT", COL.entrezgene, COl.Special.gene.List)]
##########################
TABLEs.ALL.ALL= TABLEs.ALL.ALL[,c(COl.Needs.TABLEs.ALL.ALL,"Raw.geneID","Test.XX","gene.PLOT","type.gene.PLOT", COL.entrezgene)] 
########################## 
if(length(COl.Special.gene.List)>0) for(CC.spe in COl.Special.gene.List){
	GENES.PLOT..VEC= TABLE.ALL.JOIN.plot[TABLE.ALL.JOIN.plot[,CC.spe]!="","gene.PLOT"]
	TABLEs.ALL.ALL[,CC.spe]="";  if(length(GENES.PLOT..VEC)>0) for(gene.p in GENES.PLOT..VEC) TABLEs.ALL.ALL[TABLEs.ALL.ALL[,"gene.PLOT"]==gene.p,CC.spe]="XX"
}
###########################
TABLEs.ALL.ALL.Plot= TABLEs.ALL.ALL;    DESeq2..Pvalue.NA.Values.x= c(0.5, 0.05, 0)[1]
for(DESeq.cor in c("DESeq2","DESeq2b")){ 
	TABLEs.ALL.ALL.Plot[TABLEs.ALL.ALL.Plot[,"Model"]==DESeq.cor & is.na(TABLEs.ALL.ALL.Plot[,"PValue"]),"PValue"]= DESeq2..Pvalue.NA.Values.x;  
	TABLEs.ALL.ALL.Plot[TABLEs.ALL.ALL.Plot[,"Model"]==DESeq.cor & is.na(TABLEs.ALL.ALL.Plot[,"FDR"   ]),"FDR"   ]= DESeq2..Pvalue.NA.Values.x
}

#################################################### TABLE.Models.SAVE.Com {  
####################################################
TABLE.Models.SAVE.Com= TABLE.ALL.JOIN;	  rm(TABLE.ALL.JOIN)
#####################
TABLE.Percentage.SAVE= FFF.RNAseq..JOIN.Algorithms.Percetage( TABLE.Models.SAVE.Com[TABLE.Models.SAVE.Com[,"Filters.Max"]=="",], test_RUN.Model, print.Anno)
#####################
TABLE.Models.SAVE.Com= TABLE.Models.SAVE.Com[,COL.SAVE..TABLE.ALL];     TABLE.Models.SAVE.Com= FFF.RNAseq..JOIN..Change.Columns.saving(TABLE.Models.SAVE.Com, print..RNAseq..script)

############################################################################################################ path analysis (TABLE.SAVE.TOGO)
############################################################################################################ {
print(""); print(paste("#############################################################   path analysis          ", date() )) 
Output.Kaos.or.Org.Pipe= c(1,2)[1];   FFF.RNAseq..Path.analysis..Load.MAIN.GENE.PATH.files(DIR_Annotations.Organism, TABLEs.ALL.ALL.Plot, Output.Kaos.or.Org.Pipe)

#################################################### MAster file for genelist (I use it for Monte carlo scripts)  
####################################################
SAVE.MASTER..Gene.List.per.path=c(0,1)[1];     #### stop(); #### {  SAVE.MASTER..Gene.List.per.path=c(0,1)[2]
if(SAVE.MASTER..Gene.List.per.path==1){
	######
	Go.KEGG.GeneList..MASTER= Go.KEGG.GeneLinks..MAIN.GENES..FFF;  head(Go.KEGG.GeneList..MASTER)
	######
	FILE.MASTER.X= paste(MAIN.MAIN.DIR...Inp.DataBases.ALL,"YY..RNAseq/_Gene.info/",DIR_Annotations.Organism,"/",RNAseq..Version..Gene.info,"/","MASTER..Go.KEGG.GeneList.Rdata",sep="")
	######
	save(Go.KEGG.GeneList..MASTER, file= FILE.MASTER.X )
	write.csv(Go.KEGG.GeneList..MASTER, file= gsub(".Rdata",".csv",FILE.MASTER.X,fixed=TRUE), row.names=F, na="")
	######
	PATh.UNIQUE.X= unique(Go.KEGG.GeneList..MASTER[,"PathwayID"]);   Gene.UNIQUE.X= unique(Go.KEGG.GeneList..MASTER[,"Gene"])
	print(paste(" Saving MASTERs  --> Go.KEGG.GeneLinks..MASTER  --> ",nrow(Go.KEGG.GeneList..MASTER), "           paths=  ", length(PATh.UNIQUE.X),"      genes=", length(Gene.UNIQUE.X) ))
	print(FILE.MASTER.X)
}

#####################
MAX.FDR..Res.GO..VEC= c(MAX.FDR..Res.GO, MAX.FDR..Res.GO.Extras);   MAX.FDR..Res.GO..VEC= MAX.FDR..Res.GO..VEC[order(MAX.FDR..Res.GO..VEC)]
conta.FDR=0;  
if(DIR_Annotations.Source.hg=="RSEM" | DIR_Annotations.Source.hg=="Others" | RNAseq..EMERGENCY..RUn.Pathways..Pipe==0) { 
	TABLE.SAVE.TOGO.V1= TABLEs.ALL.ALL.Plot[0,];		TABLE.SAVE.Per.TOGO= c()
	TABLE.SAVE.TOGO.V2= TABLEs.ALL.ALL.Plot[0,]
	TABLE.SAVE.TOGO.V3= TABLEs.ALL.ALL.Plot[0,]
	TABLE.SAVE.TOGO.V4= TABLEs.ALL.ALL.Plot[0,];    if(RNAseq..EMERGENCY..RUn.Pathways..Pipe==0){ print(paste("WARNING WARNING  --->   not running pathways pvalue   ----> EMERGENCY..RUn.Pathways..Pipe =", RNAseq..EMERGENCY..RUn.Pathways..Pipe ));  print("") }
} else for(MAX.FDR.GO.USE in MAX.FDR..Res.GO..VEC) { 
	###########################
	print(""); print(paste("#############################################################   path analysis         FDR=",MAX.FDR.GO.USE, "          ", date() ))   
	###########################
	CASE.Model.VEC.GO= RNAseq..Models.Alg..Save.Pipe;    CASE.Model.VEC.GO= RNAseq..Models.Alg..Max.Min.Mean;    #CASE.Model.VEC.GO= RNAseq..Models.Alg..All.Algoritms
	CASE.Model.VEC.GO= c(CASE.Model.VEC.GO, "ALL")
	########################### 
	########################### RUN.goana.Kegga
	if(MAX.FDR.GO.USE==MAX.FDR..Res.GO..Plot.QC) RUN.goana.Kegga..for.Plot= 1 else RUN.goana.Kegga..for.Plot= 0;
	
	
	
	## CASE.Model.VEC.GO= CASE.Model.VEC.GO[CASE.Model.VEC.GO!="Lima3-Voom"]
	
	##########     
	RUN.Big.Table.TOGO= 1;   PLOT..Res.GO=c(TRUE,FALSE)[1];     print..Res.GO= c(1,0)[2];    MAX.PValue.Path..Models= c(-1,0.05)[1];
	TABLE.SAVE.TOGO    = FFF.RNAseq..RUN.goana.Kegga(RUN.goana.Kegga..for.Plot, RUN.Big.Table.TOGO, TABLEs.ALL.ALL.Plot, CASE.Model.VEC.GO, test_RUN.Model, DIR_Annotations.Organism, MIN.logFC..Res.GO, MAX.PValue..Res.GO, MAX.FDR.GO.USE, MAX.PValue.Path..Models, print..Res.GO, PLOT..Res.GO) 
	TABLE.SAVE.Per.TOGO= TABLE.Percetage..TOGO..FF;  rm(TABLE.Percetage..TOGO..FF)
	### rm(TABLE.test.NEW.M1.goana.UP..FFF, TABLE.test.NEW.M1.goana.DOWN..FFF, TABLE.test.NEW.M1.Kegga.UP..FFF, TABLE.test.NEW.M1.Kegga.DOWN..FFF, TABLE.test.NEW.M2.goana..FFF, TABLE.test.NEW.M2.Kegga..FFF, GO.Kegga.TXT.filters..FFF, TXT.JOIN.filters.go.Kegga..FFF, envir = globalenv())
	### rm(GO.Kegga.Model1..FFF, GO.Kegga.Model2..FFF, GO.Kegga.TXT.filters..FFF, envir = globalenv())
	###########################
	########################### obtain genes
	if(RNAseq..Path.saving..MAX.Ngenes>0 & nrow(TABLE.SAVE.TOGO)>0) TABLE.SAVE.TOGO= TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"N"]<=RNAseq..Path.saving..MAX.Ngenes,]
	if(RNAseq..Path.saving..MIN.Ngenes>0 & nrow(TABLE.SAVE.TOGO)>0) TABLE.SAVE.TOGO= TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"N"]>=RNAseq..Path.saving..MIN.Ngenes,]
	###########################
	########################### obtain genes
	print.genes.GO=c(1,2,0)[1];     Stop_NOT.FINDING=c(1,0)[2];		  COL..RawPro.gene.GO="RawPro.gene";      COL..Gene.NAME.GO="gene.PLOT";     COL..Method.GO= "Method";  print("");
	COL.GENES.GO.KEGG= c("TNgenes.GoTo","Ngenes.GoTo","genes.GoTo","GeneID.GoTo", "N.Yes","N.Up","N.Down","N.both","genes.Up","genes.Down","genes.both"); 
	TABLE.SAVE.TOGO[,COL.GENES.GO.KEGG]= FFF.RNAseq..Obtain.genes.GO.KEGG(TABLE.SAVE.TOGO, TABLE.ALL.JOIN.plot, TABLEs.ALL.ALL.Plot, DIR_Annotations.Organism, print.genes.GO, Stop_NOT.FINDING, COL..RawPro.gene.GO, COL..Gene.NAME.GO, COL..Method.GO, MIN.logFC..Res.GO, MAX.PValue..Res.GO, MAX.FDR.GO.USE, COL.MAX.MEAN..Res.GO, CASE.Model.VEC.GO)
	ERROR.NOT.find.genes_in_path= c(ERROR.NOT.find.genes_in_path, ERROR.NOT.find.gene..FF);  rm(ERROR.NOT.find.gene..FF)	
	########################### 
	########################### add special paths and gene list
	if(length(RNAseq..Data.Special.Genes.list)>0 | length(RNAseq..Special.Pathways..List)>0) if(RNAseq..Pathanalysis..ADD.Gene.list>0) { 
		######################
		TABLE.Genes.list1= FFF.RNAseq..Special.Genes.list(RNAseq..Data.Special.Genes.list, RNAseq..Special.Pathways..List, TABLE.ALL.JOIN.plot, TABLE.SAVE.TOGO, DIR_Annotations.Organism, MIN.logFC..Res.GO, MAX.PValue..Res.GO, MAX.FDR.GO.USE, COL.MAX.MEAN..Res.GO )
		TABLE.Genes.list2= TABLE2.Special.Genes..FF;   Special.Pathways..YES= Special.Pathways..YES..FF;  rm(TABLE2.Special.Genes..FF, Special.Pathways..YES..FF)
		###################### 
		###################### Special.Pathways
		if(length(Special.Pathways..YES)>0) for(S.Path in Special.Pathways..YES) TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"RawPro.gene"]==S.Path,"User"]="Path-List"
		TEST.path.yes.ll= length(unique(TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"User"]=="Path-List","RawPro.gene"]))
		if(TEST.path.yes.ll!=length(Special.Pathways..YES)) stop(paste("MAIN ERROR --->  TEST !=length(Special.Pathways..YES) --->  ", TEST.path.yes.ll,"!=", length(Special.Pathways..YES) ))
		######################
		###################### obtain pvalues para tables
		print.pavlue=c(1,2,0)[3];   Use.Col.Max.Mean.GO= c(".Max",".Mean",".Med",".Min")[1];  MAX.FDR..GO.Genes.list= c(0.1,0.05)[1]
		TABLE.Genes.list1= FFFF..Pvalue.TABLE.Genes.list(TABLE.Genes.list1, TABLE.Models.SAVE.Com, test_RUN.Model, Use.Col.Max.Mean.GO, MIN.logFC..Res.GO, MAX.PValue..Res.GO, MAX.FDR..GO.Genes.list, print.pavlue )   ## table of special genes
		TABLE.Genes.list2= FFFF..Pvalue.TABLE.Genes.list(TABLE.Genes.list2, TABLE.Models.SAVE.Com, test_RUN.Model, Use.Col.Max.Mean.GO, MIN.logFC..Res.GO, MAX.PValue..Res.GO, MAX.FDR..GO.Genes.list, print.pavlue )   ## table of special pathways Not find in database		
		###################### 
		###################### add tables To main Data
		if(RNAseq..Pathanalysis..ADD.Gene.list==1) TABLE.SAVE.TOGO= rbind(TABLE.Genes.list1,                    TABLE.SAVE.TOGO)
		if(RNAseq..Pathanalysis..ADD.Gene.list==2) TABLE.SAVE.TOGO= rbind(                   TABLE.Genes.list2, TABLE.SAVE.TOGO)
		if(RNAseq..Pathanalysis..ADD.Gene.list==3) TABLE.SAVE.TOGO= rbind(TABLE.Genes.list1, TABLE.Genes.list2, TABLE.SAVE.TOGO)
		##########
		print("")
		print(paste(" WARNING WARNING WARNING"))
		print(paste(" WARNING WARNING WARNING --> add special paths and gene list ===> nedds to review for ensambl as --> RawPro.gene = Raw.geneID"))
		print(paste(" WARNING WARNING WARNING"))
		print("")
	}
	###########################
	###########################
	NO.Path.find.it.xx= c()
	if(length(RNAseq..Special.PathAnalalys..Terms.Name)>0) for(path.na in RNAseq..Special.PathAnalalys..Terms.Name) if(nrow(TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"Pathway"    ]==path.na,])<=0) NO.Path.find.it.xx= c(NO.Path.find.it.xx, path.na) else TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"Pathway"    ]==path.na,"User"]= "Path-Plot"
	if(length(RNAseq..Special.PathAnalalys..Terms.IDs )>0) for(path.id in RNAseq..Special.PathAnalalys..Terms.IDs ) if(nrow(TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"RawPro.gene"]==path.id,])<=0) NO.Path.find.it.xx= c(NO.Path.find.it.xx, path.id) else TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"RawPro.gene"]==path.id,"User"]= "Path-Plot"
	if(length(NO.Path.find.it.xx)>0){ NO.Path.find.it.xx= unique(NO.Path.find.it.xx);  print("");  print(paste("Warning ---> some special Pathways were not find  -->", length(NO.Path.find.it.xx)));  print(NO.Path.find.it.xx);  print("")  }
	########################### {
	########################### order Data 
	TABLE.SAVE.TOGO= TABLE.SAVE.TOGO[order(as.numeric(TABLE.SAVE.TOGO[,"pvalue.A"])),]
	##TABLE.Special.Genes= TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"User"]=="Gene-List",];       TABLE.Special.Paths= TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"User"]=="Path-List",]
	##TABLE.SAVE.TOGO= rbind(TABLE.Special.Genes, TABLE.Special.Paths, TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"User"]!="Gene-List" & TABLE.SAVE.TOGO[,"User"]!="Path-List",] ) 
	########################### 
	########################### join results 
	#print..JOIN.GO= c(1,0)[1];    #TABLE.SAVE.TOGO..ORG= TABLE.SAVE.TOGO 
	#TABLE.SAVE.TOGO = FFF.RNAseq..JOIN.Results.goana.Kegga(TABLE.SAVE.TOGO, test_RUN.Model, print..JOIN.GO);   ## TABLE.SAVE.TOGO[,COL.GENES.GO.KEGG]="";  head(TABLE.SAVE.TOGO)
	#########
	#for(col.na.xx in c("P.Down.All","P.Up.All","P.DE.All")) TABLE.SAVE.TOGO[is.na(TABLE.SAVE.TOGO[,col.na.xx]),col.na.xx]	
	########################### 
	########################### Filters column Filters
	TABLE.SAVE.TOGO[,"Filters"]="";
	if(MAX.PValue.Path..Filters1>=0 & MAX.PValue.Path..Filters1<=1) TABLE.SAVE.TOGO[as.numeric(TABLE.SAVE.TOGO[,"pvalue.A"])>MAX.PValue.Path..Filters1,"Filters"]= "LOW"
	if(MAX.PValue.Path..Filters2>=0 & MAX.PValue.Path..Filters2<=1) TABLE.SAVE.TOGO[as.numeric(TABLE.SAVE.TOGO[,"pvalue.A"])>MAX.PValue.Path..Filters2,"Filters"]= "NOT PASS"
	########################### 
	########################### ADD.GENE.LIST TO TOGO 
	print.ADD.Gene=c(1,2,0)[1]
	TABLE.SAVE.TOGO= FFF.RNAseq..ADD.GENE.LIST.TOGO(TABLE.SAVE.TOGO, TABLE.Models.SAVE.Com, RNAseq..Data.Special.Genes.list, print.ADD.Gene);   ## head(TABLE.SAVE.TOGO[,c(1:10)])
	########################### 
	########################### print special paths
	COl.print= c("RawPro.gene","Pathway","Ngenes.GoTo","##A","Filters","pvalue.A","Down.All","Up.All","DE.All","P.Down.All","P.Up.All","P.DE.All")
	TOPRINT= TABLE.SAVE.TOGO[TABLE.SAVE.TOGO[,"User"]!=0,COl.print];  if(nrow(TOPRINT)>0){ print("");  print("############ Results");  rownames(TOPRINT)=c(1:nrow(TOPRINT)); print("");  print(TOPRINT) };     print("")
	########################### 
	########################### {
	conta.FDR=conta.FDR+1
	if(conta.FDR==1) TABLE.SAVE.TOGO.V1= TABLE.SAVE.TOGO 
	if(conta.FDR==2) TABLE.SAVE.TOGO.V2= TABLE.SAVE.TOGO
	if(conta.FDR==3) TABLE.SAVE.TOGO.V3= TABLE.SAVE.TOGO
	if(conta.FDR==4) TABLE.SAVE.TOGO.V4= TABLE.SAVE.TOGO
	if(conta.FDR==5) stop()
	###########################
	TABLE.SAVE.TOGO.CSV= FFF..TABLE.TOGO..Divide.genes.Columns(TABLE.SAVE.TOGO, RNAseq..MAX.characters.saving.CSV)
	##############################
	FILE.Join.TOGO..csv..FDR= FILE.Join.TOGO..csv;   FILE.Join.TOGO..csv.Per.FDR= FILE.Join.TOGO..csv.Per
	if(MAX.FDR.GO.USE!=MAX.FDR..Res.GO..Plot.QC) {
		FILE.Join.TOGO..csv..FDR   = sub(".csv",paste("..FDR_",MAX.FDR.GO.USE,".csv",sep=""),FILE.Join.TOGO..csv    , fixed=TRUE)
		FILE.Join.TOGO..csv.Per.FDR= sub(".csv",paste("..FDR_",MAX.FDR.GO.USE,".csv",sep=""),FILE.Join.TOGO..csv.Per, fixed=TRUE)
	}
	########################### 
	XX.COLNAMES.SAVEn.TOGO= colnames(TABLE.SAVE.TOGO.CSV);    COL.REMOVE.nada.x=c();  
	######
	for(CC.tt in XX.COLNAMES.SAVEn.TOGO) if(nrow(TABLE.SAVE.TOGO.CSV[ is.na(TABLE.SAVE.TOGO.CSV[,CC.tt])                                  ,])>= nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, CC.tt)
	for(CC.tt in XX.COLNAMES.SAVEn.TOGO) if(nrow(TABLE.SAVE.TOGO.CSV[!is.na(TABLE.SAVE.TOGO.CSV[,CC.tt]) & TABLE.SAVE.TOGO.CSV[,CC.tt]=="",])>= nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, CC.tt)
	######
	STOP_ERRORS.Pathanalysis=0
	if(nrow(TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.both"]==0                               ,])==nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, "N.both","genes.both")
	if(nrow(TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"User"  ]==0                               ,])==nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, "User"  )
	if(nrow(TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.Yes" ]==TABLE.SAVE.TOGO.CSV[,"genes.A" ],])==nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, "N.Yes" ) else { 
		    TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.Yes" ]==TABLE.SAVE.TOGO.CSV[,"genes.A" ],"N.Yes" ]=NA;   STOP_ERRORS.Pathanalysis=1 }
	if(nrow(TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.Up"  ]==TABLE.SAVE.TOGO.CSV[,"Up.All"  ],])==nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, "N.Up"  ) else { 
		    TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.Up"  ]==TABLE.SAVE.TOGO.CSV[,"Up.All"  ],"N.Up"  ]=NA;   STOP_ERRORS.Pathanalysis=1 }
	if(nrow(TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.Down"]==TABLE.SAVE.TOGO.CSV[,"Down.All"],])==nrow(TABLE.SAVE.TOGO.CSV)) COL.REMOVE.nada.x=c(COL.REMOVE.nada.x, "N.Down") else { 
		    TABLE.SAVE.TOGO.CSV[TABLE.SAVE.TOGO.CSV[,"N.Down"]==TABLE.SAVE.TOGO.CSV[,"Down.All"],"N.Down"]=NA;   STOP_ERRORS.Pathanalysis=1 }			
	######
	for(CC.rm in c("Method","UP.DOWN","N1","N2", "Ont","TNgenes.GoTo",COL.REMOVE.nada.x)) XX.COLNAMES.SAVEn.TOGO= XX.COLNAMES.SAVEn.TOGO[XX.COLNAMES.SAVEn.TOGO!=CC.rm]
	TABLE.SAVE.TOGO.CSV= TABLE.SAVE.TOGO.CSV[,XX.COLNAMES.SAVEn.TOGO]
	######
	XX.COLNAMES.SAVEn.TOGO.L= gsub(paste("T",test_RUN.Model,".",sep=""),"",XX.COLNAMES.SAVEn.TOGO,fixed=TRUE);    
	XX.COLNAMES.SAVEn.TOGO.L[XX.COLNAMES.SAVEn.TOGO.L=="N"]="genes.P";		XX.COLNAMES.SAVEn.TOGO.L[XX.COLNAMES.SAVEn.TOGO.L=="Ngenes.GoTo"]="genes.D"
	colnames(TABLE.SAVE.TOGO.CSV)= XX.COLNAMES.SAVEn.TOGO.L
	###########################
	if(!file.exists(DIR_RNAseq..Results.Models)) dir.create(DIR_RNAseq..Results.Models, showWarnings = FALSE, recursive = TRUE);
	setwd(DIR_RNAseq..Results.Models);  write.csv(TABLE.SAVE.TOGO.CSV, file= FILE.Join.TOGO..csv..FDR   , row.names=F, na="");   print(paste(" Saving results --> rows=",nrow(TABLE.SAVE.TOGO.CSV  ),"      -->",FILE.Join.TOGO..csv..FDR   ))
	if(Saving..TABLE.Percentage==1 ) {  write.csv(TABLE.SAVE.Per.TOGO, file= FILE.Join.TOGO..csv.Per.FDR, row.names=T, na="");   print(paste(" Saving results --> rows=",nrow(TABLE.SAVE.Per.TOGO  ),"      -->",FILE.Join.TOGO..csv.Per.FDR))  } 
	print(DIR_RNAseq..Results.Models)
	if(STOP_ERRORS.Pathanalysis!=0) if(length(RNAseq..Path..DONT_STOP_ERROR_SUM..Test.VEC[paste(RNAseq..Path..DONT_STOP_ERROR_SUM..Test.VEC)==paste(test_RUN.Model..Pipe)])<=0) stop(" MAIN ERROR --- kaOs is stopping --> see below")
}
	
################################################################################################################ Ploting {
################################################################################################################
print("");  print(paste("########################################################################       ", Project.RNA..run, "       rows=", nrow(TABLEs.ALL.ALL.Plot), "       ", date()));         
#################
if(Run.Models..Plotting.TEST==1){ dev.off();  print("");   print(paste("  Printing plots   -->  ", FILE.PDF..TEST,"          ", date(), sep="")) }

##########################
if(Run.Models..Plotting>0) { 
	######################
	## source(paste(MAIN.MAIN.DIR...Pipeline..RNseq,"0_Function.RNAseq..DESeq.r",sep=""))
	#######
	if(!file.exists(DIR_RNAseq..Results.Models)) dir.create(DIR_RNAseq..Results.Models, showWarnings = FALSE, recursive = TRUE);
	######################
	######################
	conta.FDR=0;  #stop()
	for(MAX.FDR.GO.USE in MAX.FDR..Res.GO..VEC) {
		##################
		conta.FDR=conta.FDR+1;   TABLE.SAVE.TOGO= c()
		if(conta.FDR==1) TABLE.SAVE.TOGO= TABLE.SAVE.TOGO.V1
		if(conta.FDR==2) TABLE.SAVE.TOGO= TABLE.SAVE.TOGO.V2
		if(conta.FDR==3) TABLE.SAVE.TOGO= TABLE.SAVE.TOGO.V3
		if(conta.FDR==4) TABLE.SAVE.TOGO= TABLE.SAVE.TOGO.V4
		#######
		CASE.Model.VEC..PLOT= CASE.Model.VEC;   CASE.Model.VEC..PLOT= RNAseq..Models.Alg..Save.Pipe
		if(length(TABLE.SAVE.TOGO)<=0) { stop(paste("ERROR ERROR --> case for plotting is not run yet --> FDR =", MAX.FDR.GO.USE))
		} else {
			print("");  print(paste("################## Plotting   FDR=", MAX.FDR.GO.USE))
			###############
			###############
			Use.Col.Max.Mean.plot..VEC1= c(".Max",".Mean",".Med",".Min")[1:2]
			if(MAX.FDR.GO.USE==MAX.FDR..Res.GO..Plot.QC) for(Use.Col.Max.Mean.plot in Use.Col.Max.Mean.plot..VEC1) {
				if(RNAseq..Plooting..Normal.Cases==1) if(Use.Col.Max.Mean.plot==".Max") {
  					#########
  					if(Use.Col.Max.Mean.plot!=".Max") FILE.Join.Algorithms..PDF1.USE= sub(".pdf",paste(Use.Col.Max.Mean.plot,".pdf",sep=""), FILE.Join.Algorithms..PDF1, fixed=TRUE) else FILE.Join.Algorithms..PDF1.USE= FILE.Join.Algorithms..PDF1
  					#########
  					setwd(DIR_RNAseq..Results.Models);    pdf(FILE.Join.Algorithms..PDF1.USE);            ### cairo_pdf(FILE.Join.Algorithms..PDF1); https://www.stat.auckland.ac.nz/~paul/R/PDF/pdfEncoding.html
  					#########
  					Plot.MODEL.MAIN..VEC..USE=c()
				    FFF.RNAseq..Plotting_Models(Counts.gene.PLOT, Targets.Model, design.Model, TABLEs.ALL.ALL.Plot, TABLE.ALL.JOIN.plot, TABLE.SAVE.TOGO, CASE.Model.VEC..PLOT, Use.Col.Max.Mean.plot, RNAseq..Models..VOOM, DIR_Annotations.Organism, c(), RNAseq..Data.Special.Genes.list)  
					#########
					dev.off();  print(paste("  Printing plots   --> ", FILE.Join.Algorithms..PDF1.USE, "      ", date(), sep=""));    print(""); 
				}  
				#########
				if(length(RNAseq..Plooting..Special.Cases)>0){
					#########
  					FILE.Join.Algorithms..PDF3.USE= sub(".pdf",paste(".",Use.Col.Max.Mean.plot,".pdf",sep=""), FILE.Join.Algorithms..PDF3, fixed=TRUE)
					setwd(DIR_RNAseq..Results.Models);    pdf(FILE.Join.Algorithms..PDF3.USE); 
					#########
				    FFF.RNAseq..Plotting_Models(Counts.gene.PLOT, Targets.Model, design.Model, TABLEs.ALL.ALL.Plot, TABLE.ALL.JOIN.plot, TABLE.SAVE.TOGO, CASE.Model.VEC..PLOT, Use.Col.Max.Mean.plot, RNAseq..Models..VOOM, DIR_Annotations.Organism, RNAseq..Plooting..Special.Cases, RNAseq..Data.Special.Genes.list)
					#########
					dev.off();  print(paste("  Printing plots   --> ", FILE.Join.Algorithms..PDF3.USE, "      ", date(), sep=""));  print("");   
			}	}
			###############
			###############
			Use.Col.Max.Mean.plot..VEC2= c(".Max",".Mean",".Med",".Min")[1]
			if(PLOT.barcodeplot..GO.KEGG==1) if(nrow(TABLE.SAVE.TOGO)>0) for(Use.Col.Max.Mean.plot in Use.Col.Max.Mean.plot..VEC2){
				#########
				FILE.Join.Algorithms..PDF2..USE= FILE.Join.Algorithms..PDF2;  
				if(MAX.FDR.GO.USE!=MAX.FDR..Res.GO..Plot.QC) FILE.Join.Algorithms..PDF2..USE= sub(".pdf",paste("..FDR_",MAX.FDR.GO.USE,".pdf",sep=""),FILE.Join.Algorithms..PDF2, fixed=TRUE)
				if(length(Use.Col.Max.Mean.plot..VEC2)>1)  FILE.Join.Algorithms..PDF2..USE= sub(".pdf",paste(".",Use.Col.Max.Mean.plot,".pdf",sep=""), FILE.Join.Algorithms..PDF2..USE, fixed=TRUE)
				setwd(DIR_RNAseq..Results.Models);    pdf(FILE.Join.Algorithms..PDF2..USE);  #print("");   
				#########
				RNAseq..Plooting..barcodeplot= c("VennDiagram.GO","Path.Analysis1","barcodeplot.Gene", "barcodeplot")
				FFF.RNAseq..Plotting_Models(Counts.gene.PLOT, Targets.Model, design.Model, TABLEs.ALL.ALL.Plot, TABLE.ALL.JOIN.plot, TABLE.SAVE.TOGO, CASE.Model.VEC..PLOT, Use.Col.Max.Mean.plot, RNAseq..Models..VOOM, DIR_Annotations.Organism,RNAseq..Plooting..barcodeplot, RNAseq..Data.Special.Genes.list )
				#########
			  ##FFF.RNAseq..Plotting_Models(Counts.gene.PLOT, Targets.Model, design.Model, TABLEs.ALL.ALL.Plot, TABLE.ALL.JOIN.plot, TABLE.SAVE.TOGO, CASE.Model.VEC..PLOT, Use.Col.Max.Mean.plot, RNAseq..Models..VOOM, DIR_Annotations.Organism, "barcodeplot"     , RNAseq..Data.Special.Genes.list ) 
				#########
				dev.off();  print(paste("  Printing plots   --> ", FILE.Join.Algorithms..PDF2..USE, "      ", date(), sep=""));    #print(DIR_RNAseq..Results.Models);     print("");   
}	}	}	};	print("")
#if(exists("GO.Kegga.Model1..FFF")) rm(GO.Kegga.Model1..FFF, GO.Kegga.Model2..FFF, GO.Kegga.TXT.filters..FFF)
        
################################################################################################################ Saving CVS results  { 
################################################################################################################
print("########################################################################");    print("");  print("Saving CVS results:")

########################################################################  Counts..SAVE
########################################################################  
COL.Counts..SAVE= unique(c("##G",Extra.columns.SAVE, "RawPro.gene", "##Sam", COLUMNS.ORG..Counts.Model,"##Mean",COL.Counts.Aff,COL.Counts.UnA,"##End"))
Counts..SAVE= as.data.frame(array(NA,c(nrow(Counts.Model),length(COL.Counts..SAVE))),stringsFactors=F);      colnames(Counts..SAVE)= COL.Counts..SAVE;    rownames(Counts..SAVE)= c(1:nrow(Counts..SAVE));    
Counts..SAVE[,"RawPro.gene"]= rownames(Counts.Model);      Counts..SAVE[,grep("##",COL.Counts..SAVE,fixed=TRUE)]= "###"
#######
Counts..SAVE[,COLUMNS.ORG..Counts.Model]= Counts.Model[,COLUMNS.ORG..Counts.Model]
############ 
Run..all.extras= c(1,2)[2];  			
Counts..SAVE[,Extra.columns.SAVE]= FFF.RNAseq..Table..Annotation(Counts..SAVE, Annotation.ORG, "RawPro.gene", "RawPro.gene", Run..all.extras, Extra.columns.SAVE, print.Anno)[,Extra.columns.SAVE]
############ 
for(CASE.COL.k in c(COL.Counts.Aff,COL.Counts.UnA)) {
	if(CASE.COL.k==COL.Counts.Aff) COL.USE..kk= Targets.Model[Targets.Model[,"Test"]==RNA.SEQ.Test.VEC..Pipe[1],TXT_RNAseq..COLUMN.Samples.ID];    
	if(CASE.COL.k==COL.Counts.UnA) COL.USE..kk= Targets.Model[Targets.Model[,"Test"]==RNA.SEQ.Test.VEC..Pipe[2],TXT_RNAseq..COLUMN.Samples.ID]
	#####
	Counts..SAVE[,CASE.COL.k]=0;  
	if(length(COL.USE..kk)<=0) stop(paste(" Project can not have zero samples -->",CASE.COL.k)) 
	for(SAM.kk in COL.USE..kk)  Counts..SAVE[,CASE.COL.k]= Counts..SAVE[,CASE.COL.k] + as.numeric(Counts..SAVE[,SAM.kk])
	Counts..SAVE[,CASE.COL.k] = Counts..SAVE[,CASE.COL.k]/length(COL.USE..kk)
	#####
	Counts..SAVE[,CASE.COL.k] = as.integer(as.numeric(Counts..SAVE[,CASE.COL.k])*100)/100
}
##################### 
if(length(COl.Special.gene.List)>0){
	XX.COLNAMES.xx= colnames(Counts..SAVE);    Counts..SAVE[,COl.Special.gene.List]="";   Counts..SAVE= Counts..SAVE[,unique(c(COl.Special.gene.List,XX.COLNAMES.xx))]
	for(CC.spe in COl.Special.gene.List){
		GENES.PLOT..VEC= TABLE.Models.SAVE.Com[TABLE.Models.SAVE.Com[,CC.spe]!="","RawPro.gene"]
		Counts..SAVE[,CC.spe]="";  if(length(GENES.PLOT..VEC)>0) for(gene.p in GENES.PLOT..VEC) Counts..SAVE[Counts..SAVE[,"RawPro.gene"]==gene.p,CC.spe]="XX";   ##print(table(Counts..SAVE[,CC.spe]));  print(length(GENES.PLOT..VEC))
}	}
##################### 
Counts..SAVE= rbind(Counts..SAVE[1:2,],Counts..SAVE);     Counts..SAVE[1:2,]= "###";   Counts..SAVE[1,COLUMNS.ORG..Counts.Model]= Targets.Model[,"Test"];  Counts..SAVE[2,COLUMNS.ORG..Counts.Model]= Targets.Model[,"libType"];
Counts..SAVE[1,c(COL.Counts.Aff,COL.Counts.UnA)]= RNA.SEQ.Test.VEC..Pipe 
##################### 
ERR0R..COU= Targets.Model;  ERR0R..COU[,"TEST"]= colnames(Counts.Model);  ERR0R..COU[,"ERROR"]=0;  ERR0R..COU[ERR0R..COU[,"TEST"]!=ERR0R..COU[,TXT_RNAseq..COLUMN.Samples.ID],"ERROR"]=1;  if(nrow(ERR0R..COU[ERR0R..COU[,"ERROR"]==1,])>0) stop("MAIn ERROR target counts")

########################################################################  TABLE.Models.SAVE.Com
########################################################################     
TABLE.Models.SAVE.Com[,c(COL.Counts.Aff,COL.Counts.UnA)]= NA
NN.match= match(TABLE.Models.SAVE.Com[,"RawPro.gene"],Counts..SAVE[,"RawPro.gene"],nomatch=-1);  if(length(NN.match[NN.match<0])>0) stop("ERROR --> match RawPro.gene Counts.Aff")   ### length(NN.match);  dim(TABLE.Models.SAVE.Com);   dim(Counts..SAVE)
TABLE.Models.SAVE.Com[,c(COL.Counts.Aff,COL.Counts.UnA)]= Counts..SAVE[NN.match,c(COL.Counts.Aff,COL.Counts.UnA)]
##################### 
TABLE.Models.SAVE.Com[,c("##Mean")]= "##" 

####################################
if(!file.exists(DIR_RNAseq..Results.Models.Rdata)) dir.create(DIR_RNAseq..Results.Models.Rdata, showWarnings = FALSE, recursive = TRUE);
setwd(DIR_RNAseq..Results.Models.Rdata);     save(Counts..SAVE, TABLE.Models.SAVE.Com, TABLEs.ALL.ALL.Plot, TABLE.ALL.JOIN.plot, TABLE.Percentage.SAVE, TABLE.SAVE.TOGO, TABLE.SAVE.Per.TOGO, Desq2..vsd, dge_DESeq2, file= FILE.Join.Algorithms..Rdata);
########################
RNA.SEQ..Remove.COMLUMNS..VEC= c("Filters.Min","logFC.Min","PValue.Min","FDR.Min")
if(length(Models.remove..Pipe)>0){
	Models.remove..Pipe.id= FFF.RNAseq..JOIN..Get.Columns.Ids(Models.remove..Pipe, 0)
	RNA.SEQ..Remove.COMLUMNS..VEC= c(RNA.SEQ..Remove.COMLUMNS..VEC, paste("T",test_RUN.Model,".",Models.remove..Pipe.id,sep=""), paste("T",test_RUN.Model,".",Models.remove..Pipe.id," ",sep=""), paste("T",test_RUN.Model,".",Models.remove..Pipe.id,"  ",sep=""), paste("T",test_RUN.Model,".",Models.remove..Pipe.id,"   ",sep="") )
}
######################## 
XX.COLNAMES.SAVE  = colnames(TABLE.Models.SAVE.Com);   for(XX.rm in RNA.SEQ..Remove.COMLUMNS..VEC) XX.COLNAMES.SAVE= XX.COLNAMES.SAVE[XX.COLNAMES.SAVE!=XX.rm]
XX.COLNAMES.Counts= colnames(Counts..SAVE)
########################
for(Gene.col.x in c(RNAseq..Columns.Names..Genes.Excel,paste("T",test_RUN.Model,".D1",sep=""),paste("T",test_RUN.Model,".D1b",sep=""))){   #print(Gene.col.x) 
	#######
	if(length(XX.COLNAMES.SAVE  [XX.COLNAMES.SAVE  ==Gene.col.x])>0){ TABLE.Models.SAVE.Com[is.na(TABLE.Models.SAVE.Com[,Gene.col.x]),Gene.col.x]="";	  	TABLE.Models.SAVE.Com[TABLE.Models.SAVE.Com[,Gene.col.x]=="-",Gene.col.x]=""  } 	
	if(length(XX.COLNAMES.Counts[XX.COLNAMES.Counts==Gene.col.x])>0){ Counts..SAVE[is.na(Counts..SAVE[,Gene.col.x]),Gene.col.x]="";	  						Counts..SAVE         [Counts..SAVE         [,Gene.col.x]=="-",Gene.col.x]=""  }	
	#######
	if(length(XX.COLNAMES.SAVE  [XX.COLNAMES.SAVE  ==Gene.col.x])>0) TABLE.Models.SAVE.Com[TABLE.Models.SAVE.Com[,Gene.col.x]!="",Gene.col.x]= paste("'",TABLE.Models.SAVE.Com[TABLE.Models.SAVE.Com[,Gene.col.x]!="",Gene.col.x],"'",sep="")
	if(length(XX.COLNAMES.Counts[XX.COLNAMES.Counts==Gene.col.x])>0) Counts..SAVE         [Counts..SAVE         [,Gene.col.x]!="",Gene.col.x]= paste("'",Counts..SAVE         [Counts..SAVE         [,Gene.col.x]!="",Gene.col.x],"'",sep="")
}
######################## 
## TABLE.Models.SAVE.Com.ORG= TABLE.Models.SAVE.Com;  ### TABLE.Models.SAVE.Com= TABLE.Models.SAVE.Com.ORG 
XX.COLNAMES.SAVE2= sub(paste("T",test_RUN.Model,".",sep=""),"",XX.COLNAMES.SAVE,fixed=TRUE)
TABLE.Models.SAVE.Com= TABLE.Models.SAVE.Com[,XX.COLNAMES.SAVE];  colnames(TABLE.Models.SAVE.Com)= XX.COLNAMES.SAVE2
TABLE.Models.SAVE.Com= TABLE.Models.SAVE.Com[order(as.numeric(TABLE.Models.SAVE.Com[,"logFC.Max"]),decreasing = TRUE),]
######################## 
setwd(DIR_RNAseq..Results.Models);   print(""); 
write.csv(   Counts..SAVE                          , file= FILE.Join.Counts..csv       , row.names=F, na="");  		print(paste(" Saving results --> rows=",nrow(Counts..SAVE         ),"      -->",FILE.Join.Counts..csv       ));
write.csv(TABLE.Models.SAVE.Com[,XX.COLNAMES.SAVE2], file= FILE.Join.Algorithms.Com.csv, row.names=F, na="");  		print(paste(" Saving results --> rows=",nrow(TABLE.Models.SAVE.Com),"      -->",FILE.Join.Algorithms.Com.csv));
#########
if(Saving..TABLE.Percentage==1) {
	write.csv(TABLE.Percentage.SAVE, file= FILE.Join.Percentage..csv   , row.names=T, na="");  	print(paste(" Saving results --> rows=",nrow(TABLE.Percentage.SAVE),"      -->",FILE.Join.Percentage..csv   )); 
};  print(DIR_RNAseq..Results.Models)
#################################### testing nrow(Counts..SAVE)
if( (nrow(Counts..SAVE)-2)!=nrow(TABLE.Models.SAVE.Com)) stop(paste("MAIN ERROR -->   nrow(Counts..SAVE)-2)!=nrow(TABLE.Models.SAVE.Com) -->", (nrow(Counts..SAVE)-2),"!=",nrow(TABLE.Models.SAVE.Com) ))

############################################################################################################
############################################################################################################ testing
if(Runnning_testing_functions_inputs==1) {
	############################
	print("");   print("");   stop("Kaos is stopping --> perfect perfect --> stop test")
}

############################################################################################################
############################################################################################################
if(length(ERROR.NOT.find.genes_in_path)>0 & 66==99) {
	ERROR.NOT.find.genes_in_path= unique(ERROR.NOT.find.genes_in_path)
	print("");  print(paste("########################################################## error with GO or KEGG"))
	##############
	print(paste(" WARNING WARNING --> some paths in GO or KEGG do NOT have gene list -->", length(ERROR.NOT.find.genes_in_path)))
	if(print..RNAseq..script>0){ if(length(ERROR.NOT.find.genes_in_path)>10){ 	print("");  print(head(ERROR.NOT.find.genes_in_path, 4));    print(tail(ERROR.NOT.find.genes_in_path, 4)) } else print(ERROR.NOT.find.genes_in_path) }; print("")
	###########
	ERROR.NOT.find.genes_in_path.ALL= c(ERROR.NOT.find.genes_in_path.ALL, ERROR.NOT.find.genes_in_path)
}
############################################################################################################
############################################################################################################
print("");  print(paste("#######################################################################################################"))
            print(paste("#######################################################################################################")) 
}}}}

########################################################################################################################################################################################################################
########################################################################################################################################################################################################################
STOP.script=0
if(length(ERROR.NOT.find.genes_in_path.ALL)>0) {
	ERROR.NOT.find.genes_in_path.ALL= unique(ERROR.NOT.find.genes_in_path.ALL)
	print("");  print(paste("########################################################## error with GO or KEGG"))
	##############
	print(paste(" WARNING WARNING --> some paths in GO or KEGG do NOT have gene list -->", length(ERROR.NOT.find.genes_in_path.ALL)))
	if(print..RNAseq..script>-1){ if(length(ERROR.NOT.find.genes_in_path.ALL)>10){ 	print("");  print(head(ERROR.NOT.find.genes_in_path.ALL, 4));    print(tail(ERROR.NOT.find.genes_in_path.ALL, 4)) } else print(ERROR.NOT.find.genes_in_path.ALL) }; print("")
}
if(length(SAMPLES.QC.ERROR..ALL)>0) {
	SAMPLES.QC.ERROR..ALL= unique(SAMPLES.QC.ERROR..ALL[order(SAMPLES.QC.ERROR..ALL)]); 
	print(paste(" WARNING WARNING --> some samples dint pass QC results -->", length(SAMPLES.QC.ERROR..ALL)));  print(SAMPLES.QC.ERROR..ALL);  print("")
	STOP.script=1
}

####################################################################################################
####################################################################################################
print("");  print(paste("############################################################## END END"))
print(paste("  ****** --------- ******* --------- ******* --------- ******* --------- *******  " ))
print(paste("-----------------------------Initial time --->", Initial_TIME ))
print(paste("------------------------------ Final time --->", date()       ))
if(STOP.script==1 & Inputs.Bsub==0){ print("");  stop("kaOs is stopping") }
print("PERFECT -- Simulation complete")
print("")
####################################################################################################
####################################################################################################
####################################################################################################
setwd(WD.ORG);

####################################################################################################
####################################################################################################  edgeR (example pag 69)
############ plotMD
#### Ideally, the bulk of genes should be centred at a log-fold change of zero. This indicates that
#### any composition bias between libraries has been successfully removed. This quality check
#### should be repeated by constructing a MD plot for each sample.
####
############ plotMD
####Replicate samples from the same group cluster together in the plot, while samples from
####different groups form separate clusters. This indicates that the differences between groups
####are larger than those within groups, i.e., differential expression is greater than the variance
####and can be detected. The distance between basal samples on the left and luminal cells on
####the right is about 6 units, corresponding to a leading fold change of about 64-fold (26 = 64)
####between basal and luminal. The expression differences between virgin, pregnant and lactating
####are greater for luminal cells than for basal.
############ plotBCV
#### The square root of the common dispersion gives the coefficient of variation of biological variation
############ barcode plot
#### A barcode plot can be produced to visualize the results. Genes are ranked from left to right
#### by decreasing log-fold-change in the background of the barcode plot. Genes in the set of
#### msYgenes are represented by red bars whereas genes in the set of XiEgenes are represented by
#### blue bars. The line above the barcode shows the relative local enrichment of the vertical
#### bars in each part of the plot. This particular plot suggests that the male-specific genes tend
#### to have large positive log-fold-changes, whereas the X genes tend to have large negative
#### log-fold-changes.

####################################################################################################
####################################################################################################

############################################################ test Run.Model..edgeR with Boris
############################################################ 
if(77==99) if(length(Run.Model..edgeR)>0) if(Run.Model..edgeR>=2) {
	#########
	print(""); print(paste("#############################################################   Running edgeR 2         ", date() ))
	#########
	design2= design.Model;   #design2 = model.matrix(~ Test, data = Targets.Model);   
	######### 
    DGE_edgeR_B= DGEList(counts=Counts.Model, group=GROUP.Use)
  ##DGE_edgeR_B= DGEList(counts=Counts.Model)
    #########
  ##DGE_edgeR_B= calcNormFactors(DGE_edgeR_B, method = c("TMM","RLE","upperquartile","none")[1])  ### boris script
    DGE_edgeR_B= calcNormFactors(DGE_edgeR_B)  ### boris script  
    #########
  ##DGE_edgeR_B= estimateDisp(DGE_edgeR_B, design.Model, robust=TRUE); 
    #########
  ##DGE_edgeR_B= estimateGLMCommonDisp (DGE_edgeR_B, design2, verbose=T);  # Boris	
    DGE_edgeR_B= estimateGLMCommonDisp (DGE_edgeR_B, design2)	
    DGE_edgeR_B= estimateGLMTagwiseDisp(DGE_edgeR_B, design2)
    DGE_edgeR_B= estimateGLMTrendedDisp(DGE_edgeR_B, design2)
  	######
 	FIT1.lrt= glmFit(DGE_edgeR_B, design.Model)
 	TEST2.lrt= glmLRT(FIT1.lrt)
  	######
  ##FIT1.lrt= glmFit(DGE_edgeR_B, design2);
  ##TEST2.lrt= glmLRT(FIT1.lrt, contrast=c(0,-1))
  ##if(USE.diff.libType.DESeq>0) TEST2.lrt= glmLRT(FIT1.lrt) else TEST2.lrt= glmLRT(FIT1.lrt, contrast=c(0,-1))
    ######
	TEST2.lrt$table$FDR = p.adjust(TEST2.lrt$table$PValue, method="BH")
	#########
	FIT1.qlf= glmQLFit(DGE_edgeR_B, design.Model, robust=TRUE);
	TEST2.qlf= glmQLFTest(FIT1.qlf);
	TEST2.qlf$table$FDR = p.adjust(TEST2.qlf$table$PValue, method="BH")
	#########
	TABLE7= TEST2.lrt$table;     TABLE7.Boris1= TABLE7[TABLE7[,"PValue"]<=PRINT.PValue,];  TABLE7.Boris2= TABLE7[TABLE7[,"PValue"]<=PRINT.PValue & TABLE7[,"FDR"]<=PRINT.FDR.edgeR,];  print(paste(" RAw data=",nrow(TABLE7),"   pval<0.05=",nrow(TABLE7.Boris1),"   FDR=",nrow(TABLE7.Boris2) ));   if(print..RNAseq..script>1){ print(""); print("Results1:");  print(head(TABLE7));  print(dim(TABLE7));  print(summary(decideTests(TEST2.lrt)))  }
	TABLE8= TEST2.qlf$table;     TABLE8.Boris1= TABLE8[TABLE8[,"PValue"]<=PRINT.PValue,];  TABLE8.Boris2= TABLE8[TABLE8[,"PValue"]<=PRINT.PValue & TABLE8[,"FDR"]<=PRINT.FDR.edgeR,];  print(paste(" RAw data=",nrow(TABLE8),"   pval<0.05=",nrow(TABLE8.Boris1),"   FDR=",nrow(TABLE8.Boris2) ));   if(print..RNAseq..script>1){ print(""); print("Results2:");  print(head(TABLE8));  print(dim(TABLE8));  print(summary(decideTests(TEST2.qlf)))  }
	########
	if(nrow(TABLE7)>0) TABLE7[,"Model"]="edgeR.ltr2" else TABLE7= TABLE00
	if(nrow(TABLE8)>0) TABLE8[,"Model"]="edgeR.qlf2" else TABLE8= TABLE00
	CASE.Model.VEC= c(CASE.Model.VEC,"edgeR.ltr2","edgeR.qlf2")
	########
	dge_edgeR_B <<- DGE_edgeR_B;     fit_edgeR_B_lrt <<- FIT1.lrt;   fit_edgeR_B_qlf <<- FIT1.qlf;    res_edgeR_B_lrt <<- TEST2.lrt;    res_edgeR_B_qlf <<- TEST2.qlf
}

############################################################ edgeR pipelines (others)
############################################################ 
if(77==99) if(length(Run.Model..edgeR)>0) if(Run.Model..edgeR>=3) {
	print(paste("Running edgeR 3"))
	y3= estimateDisp(DGE); 
	y3= estimateCommonDisp(y3)
	y3= estimateTagwiseDisp(y3)
	et= exactTest(y3);  				if(print..RNAseq..script>0){ print(topTags(et));   summary(decideTests(et)) }
	topTags(et)
}

############################################################ filtering others
############################################################ 
if(77==99) if(RNAseq..QC..filtering..Counts==2){
	keep= rowSums(cpm(DGE_0)>1) >= 2
	DGE_0 = DGE_0[keep, , keep.lib.sizes=FALSE];     if(length(Annotation.Counts)>0) Annotation.Counts= Annotation.Counts[keep,]
}
########
if(77==99) if(RNAseq..QC..filtering..Counts==3){
####keep= filterByExpr(DGE_0, design.Model)
	keep= filterByExpr(DGE_0);  print(summary(keep))
	DGE_0 = DGE_0[keep, , keep.lib.sizes=FALSE];     if(length(Annotation.Counts)>0) Annotation.Counts= Annotation.Counts[keep,]
}

################################ Human-readable gene symbols can also be added to complement the Entrez identifiers for each gene (mouse)
####require(org.Mm.eg.db)
####Symbol= mapIds(org.Mm.eg.db, keys=rownames(DGE_0), keytype="ENTREZID", column="SYMBOL")
####Annotation.Counts$Human.Gene.Name= mapIds(org.Mm.eg.db, keys=rownames(DGE_0), keytype="ENTREZID", column="SYMBOL")
####table(Annotation.Counts$Human.Gene.Name==Annotation.Counts$gene.Raw)

################################################## boxplot limma
###MA <- normalizeWithinArrays(DGE_0)
###MA <- normalizeWithinArrays(DGE_0, method="none")
###boxplot(MA$M~col(MA$M),names=colnames(MA$M))
###	
###dge.NOR <- normalizeWithinArrays(dge)
###dge.NOR <- normalizeWithinArrays(logCPM)
###dge.NOR <- normalizeWithinArrays(fit1)
###cols <- Targets.Model$Test;   TEST.ID= unique(Targets.Model$Test);  if(length(TEST.ID)!=2) stop()
###cols[cols==TEST.ID[1]] <- "blue";     cols[cols==TEST.ID[2]] <- "yellow"
###boxplot(dge.NOR$M~col(dge.NOR$M),names=Targets.Model$Sample.Identifier,col=cols,xlab="boxplot",ylab="M-values")	
