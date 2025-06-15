library(tidyverse)
library(data.table)
library(colocboost)
library(argparse)
library(janitor)
library(bedr)
source('colocboost_utils.R')

###### PARSE COMMAND LINE  ARGUMENTS ######
parser <- ArgumentParser()

parser$add_argument("-vcf", "--vcf", required=TRUE, help="Path to VCF file")
parser$add_argument("-transcriptome_bed", "--transcriptome_bed", required=TRUE, help="Path to bed file")
parser$add_argument("-proteome_bed", "--proteome_bed", required=TRUE, help="Path to gwas file")
parser$add_argument("-transcriptome_covars", "--transcriptome_covars", required=TRUE, help="Path to covars file")
parser$add_argument("-proteome_covars", "--proteome_covars", required=TRUE, help="Path to gwas file")
parser$add_argument("-phenotype_id", "--phenotype_id", required=TRUE, help="phenotype_id_string")

args <- parser$parse_args()

###### SET PATHS ######## 

# Load path to bed files
transcriptome_bed_path <- args$transcriptome_bed
proteome_bed_path <- args$proteome_bed

# Load paths to covariate files
transcriptome_covars_path <- args$transcriptome_covars
proteome_covars_path <- args$proteome_covars

phenotype_id <- args$phenotype_id
VCF_file_path <- args$vcf

outfile <- paste0(phenotype_id, '_colocboost_res.RDS') 

########### LOAD FILES #########
# Load individual level data 
proteome_bed_df <- fread(proteome_bed_path)
expression_bed_df <- fread(transcriptome_bed_path)

# Print the structure of the data frames for debugging
print("Structure of proteome_bed_df:")
print(str(proteome_bed_df))

print("Structure of expression_bed_df:")
print(str(expression_bed_df))

# Read in covariate files and clean the data 
expression_covars_df <- fread(transcriptome_covars_path) %>% prep_covar_data 
proteome_covars_df <- fread(proteome_covars_path) %>% prep_covar_data 

####### PREPROCESS BED AND COVAR FILES #######

# Take expression individual data and subset to those that are in the proteomic df
print("Column names of proteome_bed_df:")
print(colnames(proteome_bed_df))

subset_transcriptomics <- expression_bed_df %>% select(any_of(colnames(proteome_bed_df)))

# Take proteomic individual data and subset to those that are in the transcriptomic df
print("Column names of expression_bed_df:")
print(colnames(expression_bed_df))

subset_proteomics <- proteome_bed_df %>% select(any_of(colnames(expression_bed_df)))

# Filter covariate data so that it corresponds to the samples with both assays 
subset_proteomic_covars <- proteome_covars_df %>% filter(row.names(.) %in% colnames(proteome_bed_df))
subset_transcriptomic_covars <- expression_covars_df %>% filter(row.names(.) %in% colnames(expression_bed_df))

########### RUN COLOCBOOST ########
# Run colocboost wrapper for protein of interest
colocboost_res <- proteome_transcriptome_coloc(phenotype_id,
                                                subset_transcriptomics,
                                                subset_proteomics,
                                                subset_transcriptomic_covars,
                                                subset_proteomic_covars,
                                                VCF_file_path
                                            )
saveRDS(colocboost_res, outfile)
