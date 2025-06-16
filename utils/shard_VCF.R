library(tidyverse)
library(data.table)
library(bedr)


tabix_query <- function(range,phenotype_id,VCF_path){
tabix_output <- tabix(range,VCF_path) 
tabix_output
}


###### PARSE ARGUMENTS ########
parser <- ArgumentParser()

parser$add_argument("-vcf", "--vcf", dest="filename", required=TRUE,
    help="Path to VCF file")
parser$add_argument("-phenotype_bed", "--phenotype_bed", dest="filename", required=TRUE,
    help="phenotype_id_string")

parser$add_argument("-padding", "--padding", dest="filename", required=TRUE,
    help="phenotype_id_string")


args <- parser$parse_args()

###### SET PATHS ######## 

# load path to bed files
phenotype_bed_path <- basename(args$phenotype_bed)
VCF_path <- basename(args$vcf)

padding <- as.numeric(args$padding)

##### PARSE DATA ######

# load in phenotyp ebed file and clean data 
# for tabix query 
phenotype_bed <- phenotype_bed_path %>% 
    fread()  %>% 
    select(1,2,3,4) %>%
    dplyr::rename('chr' = 1,'start' =2,'end' = 3,'phenotype_id' = 4)
    mutate(start = as.numeric(start) - padding,end = as.numeric(end) + padding) %>%
    mutate(start = case_when(start < 1 ~ 1,TRUE~start)) %>% 
    mutate(range  = paste0(chr,":",start,'-',end ))



