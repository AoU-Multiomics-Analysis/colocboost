library(tidyverse)
library(data.table)
library(colocboost)
library(argparse)
library(janitor)
library(bedr)
source('colocboost_utils.R')


###### PARSE COMMAND LINE  ARGUMENTS ######
parser <- ArgumentParser()

parser$add_argument("-GenotypeDosage", "--GenotypeDosage", required=TRUE, help="Path to VCF file")
parser$add_argument("-BedFile", "--BedFile", required=TRUE, help="Path to bed file")
parser$add_argument("-Covars", "--Covars", required=TRUE, help="Path to covars file")
parser$add_argument("-PhenotypeID", "--PhenotypeID", required=TRUE, help="phenotype_id_string")

args <- parser$parse_args()

###### SET PATHS ########
#Load path to bed files
BedFilePath <- args$BedFile

# Load paths to covariate files
CovarsPath <- args$Covars

PhenotypeID <- args$PhenotypeID

outfile <- paste0(PhenotypeID, '_colocboost_res.RDS') 


####### LOAD FILES ##########

# load bed file and prep covariate data 
BedData <- fread(BedFilePath)
CovarsData <- CovarsPath %>% prep_covar_data 



