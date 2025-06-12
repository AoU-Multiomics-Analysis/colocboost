library(tidyverse)
library(data.table)
library(colocboost)

source('colocboost_utils.R')

###### PARSE COMMAND LINE  ARGUMENTS ######

parser$add_argument("-vcf", "--vcf", dest="filename", required=TRUE,
    help="Path to VCF file")
parser$add_argument("-transcriptome_bed", "--transcriptome_bed", dest="filename", required=TRUE,
    help="Path to bed file")
parser$add_argument("-proteome_bed", "--proteome_bed", dest="filename", required=TRUE,
    help="Path to gwas file")
parser$add_argument("-transcriptome_covars", "--transcriptome_covars", dest="filename", required=TRUE,
    help="Path to covars file")
parser$add_argument("-proteome_covars", "--proteome_covars", dest="filename", required=TRUE,
    help="Path to gwas file")
parser$add_argument("-phenotype_id", "--phenotype_id", dest="filename", required=TRUE,
    help="phenotype_id_string")


args <- parser$parse_args()

###### LOAD PATHS ######## 

# load path to bed files
transcriptome_bed_path <- basename(args$transcriptome_bed)
proteome_bed_path <- basename(args$proteome_bed)

# load paths to covariate files
transcriptome_covars_path <- basename(args$transcriptome_covars)
proteome_covars_path <- basename(args$proteome_covars)


phenotype_id <- args$phenotype_id
VCF_file_path <- basename(args$vcf)



##### LOAD FILES #########
proteome_bed_df <- fread(proteome_bed_path)
expression_bed_df <- fread(transcriptome_bed_path)

expression_covars_df <- fread(transcriptome_covars_path) %>% prep_covar_data 
proteome_covars_df <- fread(proteome_covars_path) %>% prep_covar_data 

