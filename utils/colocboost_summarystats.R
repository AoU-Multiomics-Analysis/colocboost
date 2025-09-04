message('Begin colocboost summary stats')
library(tidyverse)
library(data.table)
library(colocboost)
library(argparse)
library(janitor)
library(bedr)
source('/tmp/colocboost_utils.R')


###### PARSE COMMAND LINE  ARGUMENTS ######
parser <- ArgumentParser()

parser$add_argument("-GenotypeDosage", "--GenotypeDosage", required=TRUE, help="Path to VCF file")
parser$add_argument("-BedFile", "--BedFile", required=TRUE, help="Path to bed file")
parser$add_argument("-Covars", "--Covars", required=TRUE, help="Path to covars file")
parser$add_argument("-PhenotypeID", "--PhenotypeID", required=TRUE, help="phenotype_id_string")
parser$add_argument("-SumstatsGWAS", "--SumstatsGWAS", required=TRUE, nargs = "+")



args <- parser$parse_args()


clean_GWAS_data <- function(input_GWAS_data) {
GWAS_name <- str_remove(basename(input_GWAS_data),'_munged_summary_statistics.tsv.gz')

GWAS_dat <- fread(input_GWAS_data) %>% 
                dplyr::rename('chromosome' = 2,
                                'base_pair_location' = 3,
                                'beta' =6,
                                'n' = 'N',
                                'standard_error' = 'SE')  %>%
                mutate(chromosome = paste0('chr',chromosome))

setNames(list(GWAS_dat),GWAS_name)
}

###### SET PATHS ########
message('Loading command line arguments')

#Load path to bed files
BedFilePath <- args$BedFile

# Load paths to covariate files
CovarsPath <- args$Covars

PhenotypeID <- args$PhenotypeID

Outfile <- paste0(PhenotypeID, '_colocboost_res.RDS') 

DosageFile <- args$GenotypeDosage

SumstatsGWAS <- args$SumstatsGWAS

####### LOAD FILES ##########

# load bed file and prep covariate data 
message('Loading Bedfile')
BedData <- read_table(BedFilePath)

message('Loading covariates')
CovarsData <- CovarsPath %>% prep_covar_data 

message('Preprocessing colocboost data')
ColocboostObj <- preprocess_gene_coloc_boost(PhenotypeID,BedData,DosageFile,CovarsData)

message('Loading GWAS data')
ListGWAS <- SumstatsGWAS %>%
  map(clean_GWAS_data) %>%
  flatten()

message('Extracting QTL variants from GWAS')
SumstatData <- ColocboostObj %>%
    extract_GWAS_data_list(.,ListGWAS)


message('Colocboost begin')
ColocboostResult<- colocboost(X = ColocboostObj$resid_genotype_matrix,
           Y = ColocboostObj$resid_phenotype_vec,
           LD = ColocboostObj$LD_matrix,
           sumstat = SumstatData
           )
message('Saving colocboost results')
saveRDS(ColocboostResult,Outfile)

