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

###### SET PATHS ######## 

# load path to bed files
transcriptome_bed_path <- basename(args$transcriptome_bed)
proteome_bed_path <- basename(args$proteome_bed)

# load paths to covariate files
transcriptome_covars_path <- basename(args$transcriptome_covars)
proteome_covars_path <- basename(args$proteome_covars)


phenotype_id <- args$phenotype_id
VCF_file_path <- basename(args$vcf)

outfile <- paste0(phenotype_id,'_colocboost_res.RDS') 
########### LOAD FILES #########
# load individual level data 
proteome_bed_df <- fread(proteome_bed_path)
expression_bed_df <- fread(transcriptome_bed_path)


# read in covariate files and clean the data 
expression_covars_df <- fread(transcriptome_covars_path) %>% prep_covar_data 
proteome_covars_df <- fread(proteome_covars_path) %>% prep_covar_data 


####### PREPROCESS BED AND COVAR FILES #######

# take expression individual data and subset to those that  are in the proteomic df
subset_transcriptomics <- expression_bed_df %>% select(any_of(colnames(proteomic_bed_df)))

# take proteomic individual data and subset to those that  are in the transcriptomic df
subset_proteomics <- proteomic_bed_df %>% select(any_of(colnames(expression_bed_df)))

# filter covariate data so that it corresponds to the samples with both assays 
subset_proteomic_covars <- proteomic_covars_df %>% filter(row.names(.) %in% colnames(proteomic_bed_df))
subset_transcriptomic_covars <- expression_covars_df %>% filter(row.names(.) %in% colnames(expression_bed_df))


########### RUN COLOCBOOST ########
# run colocboost wrapper for protein of interest
colocboost_res <- proteome_transcriptome_coloc(phenotype_id,
                                                subset_transcriptomics,
                                                subset_proteomics,
                                                subset_transcriptomic_covars,
                                                subset_proteomic_covars,
                                                VCF_path
                                            )
saveRDS(colocboost_res,outfile)


