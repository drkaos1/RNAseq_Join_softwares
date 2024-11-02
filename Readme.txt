
##########################################################################################################
##########################################################################################################
##########################################################################################################
RNA-seq-Join_softwares
Code to integrate 3 RNAseq softwares

This code integrate the results for 3 RNA-seq softwares (limma, edgeR and DESeq2) and these are the steps to run this code

##########################################################################################################
########################################################################################################## run example
##########################################################################################################

This code is ready to use our example "Example..P1". 

a) Unzip the file "RNA_seq_Join_softwares.zip" in one directory (working director). 
b) Change your working director in "RR3.Integrete.RNAseq.Soft.Algorithms.r (input= MAIN_WROKING_Directory)
c) Run the code in R. Make sure you have installed the all the libraries.


##########################################################################################################
########################################################################################################## run code with your project
##########################################################################################################

############################################################ Step 1
############################################################
1.1) Unzip the file "RNA_seq_Join_softwares.zip" in one directory (working director). 
1.2) You need to install the following libraries, see below for versions

library(limma);  
library(edgeR);   
library(Rsubread);	
library(biomaRt) 
library("DESeq");   
library("DESeq2");  
library("RColorBrewer");  
library("gplots");   
library("pheatmap") 

##################### Versions (sessionInfo)
#####################

sessionInfo() 

R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /n/groups/markianos/klaus/88_Programs/R_Versions/R-3.6.3/lib/libRblas.so
LAPACK: /n/groups/markianos/klaus/88_Programs/R_Versions/R-3.6.3/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pheatmap_1.0.12             gplots_3.1.0                RColorBrewer_1.1-2          DESeq2_1.26.0               SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
 [7] BiocParallel_1.20.1         matrixStats_0.57.0          GenomicRanges_1.38.0        GenomeInfoDb_1.22.1         DESeq_1.38.0                lattice_0.20-41            
[13] locfit_1.5-9.4              org.Mm.eg.db_3.10.0         GO.db_3.10.0                AnnotationDbi_1.48.0        IRanges_2.20.2              S4Vectors_0.24.4           
[19] Biobase_2.46.0              BiocGenerics_0.32.0         biomaRt_2.45.9              Rsubread_2.0.1              edgeR_3.28.1                limma_3.42.2               

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1       ellipsis_0.3.2         htmlTable_2.1.0        XVector_0.26.0         base64enc_0.1-3        rstudioapi_0.11        farver_2.0.3          
 [8] bit64_4.0.5            mvtnorm_1.1-1          apeglm_1.8.0           xml2_1.3.2             splines_3.6.3          geneplotter_1.64.0     knitr_1.30            
[15] Formula_1.2-3          annotate_1.64.0        cluster_2.1.0          dbplyr_1.4.4           png_0.1-7              compiler_3.6.3         httr_1.4.3            
[22] backports_1.1.10       assertthat_0.2.1       Matrix_1.2-18          htmltools_0.5.0        prettyunits_1.1.1      tools_3.6.3            coda_0.19-4           
[29] gtable_0.3.0           glue_1.4.2             GenomeInfoDbData_1.2.2 dplyr_1.0.2            rappdirs_0.3.1         Rcpp_1.0.5             bbmle_1.0.23.1        
[36] vctrs_0.3.4            xfun_0.18              stringr_1.4.0          lifecycle_0.2.0        gtools_3.8.2           statmod_1.4.34         XML_3.99-0.3          
[43] zlibbioc_1.32.0        MASS_7.3-53            scales_1.1.1           hms_0.5.3              curl_4.3.2             memoise_1.1.0          gridExtra_2.3         
[50] ggplot2_3.3.2          emdbook_1.3.12         bdsmatrix_1.3-4        rpart_4.1-15           latticeExtra_0.6-29    stringi_1.5.3          RSQLite_2.2.1         
[57] genefilter_1.68.0      checkmate_2.0.0        caTools_1.18.0         rlang_0.4.8            pkgconfig_2.0.3        bitops_1.0-6           purrr_0.3.4           
[64] labeling_0.3           htmlwidgets_1.5.2      bit_4.0.4              tidyselect_1.1.0       plyr_1.8.6             magrittr_1.5           R6_2.4.1              
[71] generics_0.0.2         Hmisc_4.4-1            DBI_1.1.0              pillar_1.4.6           foreign_0.8-75         survival_3.2-7         RCurl_1.98-1.2        
[78] nnet_7.3-14            tibble_3.0.4           crayon_1.3.4           KernSmooth_2.23-17     BiocFileCache_1.10.2   jpeg_0.1-8.1           progress_1.2.2        
[85] grid_3.6.3             data.table_1.13.0      blob_1.2.1             digest_0.6.25          xtable_1.8-4           numDeriv_2016.8-1.1    openssl_1.4.3         
[92] munsell_0.5.0          askpass_1.1 


############################################################ Step 2
############################################################
2.1) Define the name of you project, example "Example..P1"
2.2) Define your project in "RR1.RNAseq..Definitions.r" (search for Example..P1)
2.3) Define your test    in "RR1.RNAseq..Tests.r"       (search for Example..P1)

Please see example "Example..P1" for more details

############################################################ Step 3
############################################################
3) Your project needs 2 inputs. They need to be located inside your working director
	a) Pedigree
	b) fastq process data with annotations in R
	
Please see example "Example..P1" for more details and location

############################################################ Step 4
############################################################
4) Run the code using "RR3.Integrete.RNAseq.Soft.Algorithms.r"
	4.1) Open R3.Integrete.RNAseq.Soft.Algorithms.r
	4.2) Add your project name in "Project.RNA..run"
	4.3) Select your Annotations.Source ("Ensembl:Ch37","NCBI:Ch37","UCSC:hg19","UCSC:hg38")
	4.4) Change your working director in "MAIN_WROKING_Directory"
	4.5) Run your code in R
	














 





