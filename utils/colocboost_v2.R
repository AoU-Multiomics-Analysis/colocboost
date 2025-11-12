library(pecotmr)
library(colocboost)
library(optparse)
library(tidyverse)
library(data.table)

########### WRAPPER FUNCTIONS ##############
# Takes in the phenotype table and paths to bedfile and covariates 
# and parses arguments for load_multitask_regional_data.
# Phenotype taable should contain the columns 
# Region, gene_id, Expression, Splicing, Protein 
# where the last three columns correspond to the labels of the 
# gene in each of the three bed files
extract_parameters <- function(phenotype_table ,CovarsDf,BedFileDf) {
    long_table <- phenotype_table %>% 
        pivot_longer(!c(gene_id,region)) %>%
        arrange(name) %>% 
        filter(!is.na(value)) %>% 
        left_join(BedFileDf,by = c('name' = 'Assay')) %>%
        left_join(CovarsDf,by = c('name' = 'Assay'))
    
    region <- long_table %>% distinct(region) %>% pull(region)
    covar_paths <- long_table %>% pull(CovarsPath)
    bed_paths <- long_table %>% pull(BedPath)
    phenotype_names <- long_table %>% pull(value)
    modalities <- long_table %>% pull(name)
    match_geno_pheno <- c(rep(1,length(phenotype_names)))
    
    output <- list(
                   Region = region,
                   RegionNames = phenotype_names,
                   ConditionList = modalities,
                   MatchGenoPheno = match_geno_pheno,
                   PhenotypeList = bed_paths,
                   CovarsList = covar_paths
                  )
    output
}

extract_regional_data <- function(RegionalDataParameters,GenotypeFile) {
    region_name_col <- 4
    keep_indel = TRUE

    maf_cutoff = 0
    mac_cutoff = 0 
    xvar_cutoff = 0
    imiss_cutoff = 1

    match_geno_pheno <- RegionalDataParameters$MatchGenoPheno
    conditions_list_individual <- RegionalDataParameters$ConditionList
    phenotype_list <- RegionalDataParameters$PhenotypeList
    covariate_list <- RegionalDataParameters$CovarsList
    region <- RegionalDataParameters$Region
    association_window <- region
    extract_region_name <- lapply(RegionalDataParameters$RegionNames ,function(i) c(i))

    region_data_individual <- load_multitask_regional_data(
        region = region,
        genotype_list = GenotypeFile,
        phenotype_list = phenotype_list,
        covariate_list = covariate_list,
        conditions_list_individual = conditions_list_individual,
        match_geno_pheno = match_geno_pheno,
        association_window = association_window,
        region_name_col = region_name_col,
        extract_region_name = extract_region_name,
        keep_indel = keep_indel,
        keep_samples = SampleList,
        maf_cutoff = maf_cutoff,
        mac_cutoff = mac_cutoff,
        xvar_cutoff = xvar_cutoff,
        imiss_cutoff = imiss_cutoff
    )
    region_data_individual
}




######## LOAD COMMAND LINE ARGUMENTS ########
option_list <- list(
  optparse::make_option(c("--PhenotypeTable"), type="character", default=NULL,
    help="Table of phenotypes across modalities and windows to be tested", metavar="type"),
  optparse::make_option(c("--ExpressionMatrix"), type="character", default=NULL,
    help="Bed file containing expression values", metavar="type"),
  optparse::make_option(c("--SplicingMatrix"), type="character", default=NULL,
    help="Bed file containing splicing values", metavar="type"),
  optparse::make_option(c("--ProteinMatrix"), type="character", default=NULL,
    help="Bed file containing protein values", metavar="type"),
  optparse::make_option(c("--ExpressionCovars"), type="character", default=NULL,
    help="covars for expression matrix", metavar="type"),
  optparse::make_option(c("--SplicingCovars"), type="character", default=NULL,
    help="covars for splicing matrix", metavar="type"),
  optparse::make_option(c("--ProteinCovars"), type="character", default=NULL,
    help="covars for protein matrix", metavar="type"),
  optparse::make_option(c("--PlinkBedGenotypes"), type="character", default=NULL,
    help="path to plink genotypes", metavar="type"),
  optparse::make_option(c("--SampleList"), type="character", default=NULL,
    help="path to list of samples to use", metavar="type")
  )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))



###### PARSE COMMAND LINE ARGUMENTS ##########

PhenotypeTable <- fread(opt$PhenotypeTable)

ExpressionCovarsPath <- opt$ExpressionCovars
SplicingCovarsPath <- opt$SplicingCovars
ProteinCovarsPath <- opt$ProteinCovars
CovarsPaths <- c(ExpressionCovarsPath,SplicingCovarsPath,ProteinCovarsPath)


ExpressionMatrixPath <- opt$ExpressionMatrix
SplicingMatrixPath <- opt$SplicingMatrix
ProteinMatrixPath <- opt$ProteinMatrix
PhenotypePaths <- c(ExpressionMatrixPath,SplicingMatrixPath,ProteinMatrixPath)
ModalitiesOrder <- c('Expression','Splicing','Protein')


BedFileDf <- data.frame(BedPath=PhenotypePaths,Assay = ModalitiesOrder )
CovarsDf <- data.frame(CovarsPath=CovarsPaths,Assay = ModalitiesOrder )


SampleList <- fread(opt$SampleList) %>% dplyr::rename('value' = 1) %>% pull(value)
NumberSamples <- SampleList %>% length()
message(paste0('Number of samples included:', NumberSamples))

GenotypeFile <- sub("\\.[^.]*$", "", opt$PlinkBedGenotypes)
message(paste0('Using ',GenotypeFile,' as plink input'))

message('Extracting regional parameters')
Parameters <- PhenotypeTable %>% 
                       extract_parameters(CovarsDf,BedFileDf)

message('Extracting regional data')
RegionalData <- Parameters %>% 
                        extract_regional_data(GenotypeFile)

maf_cutoff = 0
pip_cutoff_to_skip_ind = rep(0, length(Parameters$PhenotypeList))
qc_method = "rss_qc" 


colocboost_results <- colocboost_analysis_pipeline(
    RegionalData,
    maf_cutoff = maf_cutoff, 
    pip_cutoff_to_skip_ind = pip_cutoff_to_skip_ind,
    qc_method = qc_method,
    xqtl_coloc = TRUE,
    joint_gwas = TRUE,
    separate_gwas = TRUE,
    overlap = FALSE,
    #check_null = .001,
    check_null_max_ucos = 0.0015
)

saveRDS(colocboost_results,'test.rds')




############ SET PATHS ######################
#covariate_list <- c('/home/jupyter/RNA_proteome_colocboost/covars/EUR_QTL_covariates.tsv',
		 #'/home/jupyter/RNA_proteome_colocboost/covars/EUR_splicing_QTL_covariates.tsv'
		 #)

#phenotype_list <- c('/home/jupyter/RNA_proteome_colocboost/phenotype/EUR.expression.bed.gz',
			#'/home/jupyter/RNA_proteome_colocboost/phenotype/EUR.splicing.bed.gz'
		   #)
#region <- 'chr1:206645688-208645689'
#association_window <- region
#extract_region_name <- list(c('ENSG00000307252.1'),
				#c('chr1:207646968:207647317:clu_4655_+:ENSG00000307252.1')
			   #)
#conditions_list_individual <- c('Expression','Splicing')
#match_geno_pheno <- c(1,1)

#genotype_list <- '/home/jupyter/RNA_proteome_colocboost/genotype/expression/EUR'
#region_name_col <- 4
#keep_indel = TRUE

#maf_cutoff = 0
#mac_cutoff = 0 
#xvar_cutoff = 0
#imiss_cutoff = 1

####### EXTRACT REGIONAL DATA ###############


# More advanced parameters see pecotmr::load_multitask_regional_data()
#region_data_individual <- load_multitask_regional_data(
    #region = region,
    #genotype_list = genotype_list,
    #phenotype_list = phenotype_list,
    #covariate_list = covariate_list,
    #conditions_list_individual = conditions_list_individual,
    #match_geno_pheno = match_geno_pheno,
    #association_window = association_window,
    #region_name_col = region_name_col,
    #extract_region_name = extract_region_name,
    #keep_indel = keep_indel,
    #keep_samples = SampleList,
    #maf_cutoff = maf_cutoff,
    #mac_cutoff = mac_cutoff,
    #xvar_cutoff = xvar_cutoff,
    #imiss_cutoff = imiss_cutoff
#)
#saveRDS(region_data_individual,'test_region_data.rds')
#region_data_individual <- readRDS('test_region_data.rds')

#maf_cutoff = 0
#pip_cutoff_to_skip_ind = rep(0, length(phenotype_list))
#qc_method = "rss_qc" 


#colocboost_results <- colocboost_analysis_pipeline(
    #region_data_individual,
    #maf_cutoff = maf_cutoff, 
    #pip_cutoff_to_skip_ind = pip_cutoff_to_skip_ind,
    #qc_method = qc_method,
    #xqtl_coloc = TRUE,
    #joint_gwas = TRUE,
    #separate_gwas = TRUE,
    #min_cluster_corr = 0.95,
    #overlap = FALSE,
    ##check_null = .001,
    #check_null_max_ucos = 0.0015
#)

#summary <- get_colocboost_summary(colocboost_results):w
#saveRDS(colocboost_results,'test.rds')


