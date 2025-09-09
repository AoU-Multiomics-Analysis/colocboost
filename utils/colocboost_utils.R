### Preprocessing functions
# helper function to load munged sumstats 
# and format for MR 
load_gwas_data <- function(GWAS_path,bed_range){
#GWAS_dat <- fread(GWAS_path)
#GWAS_dat_cols <- colnames(GWAS_dat)

GWAS_dat_cols <- strsplit(readLines(GWAS_path,n =1 ),split ='\t') %>% unlist()

GWAS_dat <- tabix(bed_range,GWAS_path)
colnames(tabix_res) <- GWAS_dat_cols


message(paste0('GWAS columns: ',GWAS_dat_cols,'\n'))

if (!'FRQ' %in% GWAS_dat_cols){
message('Missing FRQ column in GWAS data\n Adding column in as NA')
GWAS_dat$FRQ <- NA

}



if (!'SE' %in% GWAS_dat_cols & 'OR' %in% GWAS_dat_cols){
message('SE measurement is missing, computing from OR and P value')
GWAS_dat$SE <- get_se(GWAS_dat$OR,GWAS_dat$P)
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'OR')

    
} else if (!'SE' %in% GWAS_dat_cols & 'BETA' %in% GWAS_dat_cols){
messasge('SE measurement is missing, computing from BETA and P value')
GWAS_dat$SE <- TwoSampleMR::get_se(GWAS_dat$BETA,GWAS_dat$P)
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'BETA')

} else if ('SE' %in% GWAS_dat_cols & 'BETA' %in% GWAS_dat_cols){
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'BETA')
}else if ('SE' %in% GWAS_dat_cols & 'OR' %in% GWAS_dat_cols){
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'OR')
}

GWAS_dat
}


prep_covar_data <- function(QTL_covar_input){
covar_df <- QTL_covar_input %>%
    janitor::row_to_names(1) %>%
    filter(!is.na(ID)) %>%
    filter(!str_detect(ID,'batch')) %>%
    t() %>%
    data.frame() %>% rownames_to_column('ID') %>%
    janitor::row_to_names(1) %>% data.frame()

rownames(covar_df) <- NULL
output <- covar_df  %>%
    column_to_rownames('ID')
output
}

partialize_data <- function(covars_df,gt_matrix,pheno_vec) {

covariates_matrix <- covars_df %>% data.matrix()
hat = diag(nrow(covariates_matrix)) - covariates_matrix %*% solve(crossprod(covariates_matrix)) %*% t(covariates_matrix)
expression_vector = hat %*% phenotype_vec
subset_gt_matrix = data.matrix(gt_matrix[rownames(pheno_vec),])
gt_hat = hat %*% subset_gt_matrix 
output <- list(pheno_vec = expression_vector,gt_hat = gt_hat)
output 
    
}


# residualizes phenotype data on covariates
residualize_molecular_phenotype_data <- function(phenotype_vector,covars_df){
joint_data <- phenotype_vector %>%
        data.frame() %>%
        rownames_to_column('ID') %>%
        left_join(covars_df %>% rownames_to_column('ID'),by = 'ID') %>%
        column_to_rownames('ID') %>%
        mutate(across(everything(),~as.numeric(.)))
regression_model <- lm(value ~ .,data = joint_data)
residual_data <- bind_cols(ID = rownames(joint_data),value = resid(regression_model)) %>% column_to_rownames('ID')
residual_data
}


residualize_genotype <- function(outcome_vector,predictor_df){
joint_data <- bind_cols(outcome = outcome_vector,predictor_df) %>% mutate(across(everything(),~as.numeric(.)))
regression_model <- lm(outcome_vector ~ .,data = joint_data)
residual_data <- resid(regression_model)
residual_data

}


# pulls out expression values for a given gene from a bed file
extract_phenotype_vector <- function(phenotype_bed_file,protein_name){
phenotype_data <- phenotype_bed_file %>% filter(gene_id == protein_name) %>%
        select(-`#chr`,-start,-end) %>%
        pivot_longer(!gene_id) %>%
        column_to_rownames('name') %>%
        select(-gene_id) %>%
        data.matrix()
phenotype_data
}

# uses the phenotype bed file to extract cis-window for gene of interest and
# gets all variants within that window from VCF file using tabix
extract_genotype_vector <- function(phenotype_bed_file,protein_name,tabix_path){
require(bedr)
    
tabix_header <- strsplit(readLines(tabix_path,n =1 ),split ='\t') %>% unlist()

input_range <- phenotype_bed_file %>%
    filter(gene_id == protein_name) %>%
    mutate(start = start -1000000,end = end + 1000000) %>%
    mutate(start = case_when(start  < 1 ~ 1,TRUE ~ start) ) %>%
    mutate(range = paste0(`#chr`,':',start,'-',end)) %>%
    pull(range)

# Print debug information
print(paste("Input range:", input_range))
print(paste("Tabix path:", tabix_path))
    
tabix_res <- tabix(input_range[1],tabix_path)
colnames(tabix_res) <- tabix_header
tabix_res
}


clean_genotype_data_dosage <- function(tabix_output){
message('Cleaning genotype dosage data')
variant_metadata <- tabix_output %>%
        select(CHROM,POS,REF,ALT) %>%
        mutate(ID = paste0(CHROM,"-",POS,'-',REF,'-',ALT)) %>% select(ID)

genotype_matrix <- tabix_output %>%
        mutate(ID = paste0(CHROM,"-",POS,'-',REF,'-',ALT)) %>% 
        select(-CHROM,-POS,-ID,-REF,-ALT) %>%
        mutate(across(everything(),~as.integer(.))) %>% 
        zoo::na.aggregate(fun = mean)

output_data <- bind_cols(variant_metadata,genotype_matrix) %>%
        column_to_rownames('ID') %>%
        t() %>%
        data.frame() %>%
        mutate(across(everything(),~scale(.,center = TRUE,scale = FALSE))) %>%
        dplyr::rename_with(~str_replace_all(.,'\\.','-'))
output_data
}


clean_genotype_data <- function(tabix_output){
variant_metadata <- tabix_output %>%
        select(CHROM,POS,REF,ALT) %>%
        mutate(ID = paste0(CHROM,"-",POS,'-',REF,'-',ALT)) %>% select(ID)
genotype_matrix <- tabix_output %>%
        select(-CHROM,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT) %>%
        mutate(across(everything(),~str_remove(.,':.*')))  %>%
        mutate(across(everything(),~case_when(. == '0/0' ~ 0, . == '1/0' ~1,. == '0/1' ~1,. == '1/1' ~ 2,
                                              . == '0|0' ~ 0, . == '1|0' ~1,. == '0|1' ~1,. == '1|1' ~ 2))) %>% 
        zoo::na.aggregate()

output_data <- bind_cols(variant_metadata,genotype_matrix) %>%
        column_to_rownames('ID') %>%
        t() %>%
        data.frame() %>%
        mutate(across(everything(),~scale(.))) %>%
        dplyr::rename_with(~str_replace_all(.,'\\.','-'))
output_data
}

# extracts variant metadata from genotype matrix to be used
# to query GWAS data
extract_variant_metadata <- function(genotype_matrix){
variant_data <- data.frame(variants = colnames(genotype_matrix)) %>%
        separate(variants, into = c('chromosome','pos','ref','alt')) %>%
        mutate(pos = as.numeric(pos))
variant_data

}

extract_GWAS_data <- function(variant_metadata,GWAS_sumstats){
GWAS_data <- variant_metadata %>%
    mutate(pos = as.numeric(pos)) %>%
    mutate(variant = paste0(chromosome,'-',pos,'-',ref,'-',alt)) %>%
    left_join(GWAS_sumstats,by = c('chromosome','pos' = 'base_pair_location') ) %>%
    select(variant,beta,standard_error,n) %>%
    dplyr::rename('sebeta' = 'standard_error')
GWAS_data
}



# preps individual level data for coloc boost
preprocess_gene_coloc_boost <- function(phenotype_id,phenotype_bed,VCF_path,covars_df){

message(phenotype_id)

message('Extracting genotype data')
genotype_data  <- phenotype_bed %>%
    extract_genotype_vector(phenotype_id,VCF_path)
genotype_matrix <- genotype_data %>% clean_genotype_data_dosage
#filtered_genotype_matrix <- genotype_matrix[ , colSums(is.na(genotype_matrix)) == 0]
#filtered_genotype_matrix <- genotype_matrix[rowSums(is.na(genotype_matrix)) == 0 , ]



message('Extracting phenotype data')
phenotype_vec <- phenotype_bed %>% extract_phenotype_vector(phenotype_id)

residualized_phenotype_vec <- phenotype_vec %>%  residualize_molecular_phenotype_data(covars_df)
variant_metadata <-  extract_variant_metadata(genotype_matrix)

phenotype_range <- phenotype_bed %>% 
        filter(gene_id == phenotype_id) %>%
        dplyr::select(`#chr`,start,end) %>%
        transmute(range = paste0(`#chr`,':',start,'-',end)) %>% 
        pull(range)
 

message('Computing LD')
LD_matrix <- get_cormat(genotype_matrix)


message('Creating output object')
output <- list(LD_matrix = LD_matrix,
               variant_metadata = variant_metadata,
               resid_phenotype_vec =residualized_phenotype_vec,
               phenotype_vec = phenotype_vec,
               genotype_matrix = genotype_matrix,
               genotype_data = genotype_data,
               cis_window = phenotype_range,
               name = phenotype_id)
output
}

# takes list of summary stats of runs extract_GWAS_data function to get all
# summary stats of interest
extract_GWAS_data_list <- function(preprocessed_colocboost,sumstats_list){
variant_metadata <- preprocessed_colocboost$variant_metadata
GWAS_out <- sumstats_list %>%
                map(~extract_GWAS_data(variant_metadata,.) %>% na.omit())
names(GWAS_out) <- names(sumstats_list)
GWAS_out
}


wrap_colocboost_list <- function(colocboost_preproc_object,sumstats_list){
sumstat_data <- colocboost_preproc_object %>%
    extract_GWAS_data_list(.,GWAS_summary_stats)
try(res <- colocboost(X = colocboost_preproc_object$genotype_matrix,
           Y = colocboost_preproc_object$resid_phenotype_vec,
           LD = colocboost_preproc_object$LD_matrix,
           sumstat = sumstat_data
           ))
out <- colocboost_preproc_object
out$res <- res
try(out$summary <- res %>% get_cos_summary())
out$traits <- names(sumstat_data)
out
}


proteome_transcriptome_coloc <- function(phenotype_id,
                                         transcriptome_bed,
                                         proteome_bed,
                                         transcriptome_covar,
                                         proteome_covar,
                                         VCF_path){

message(paste0('Protein/Gene running:',phenotype_id))
split_phenotype <- unlist(str_split(phenotype_id,'_'))[2]
transcriptomic_dat <- extract_phenotype_vector(transcriptome_bed,split_phenotype)
proteomic_dat <- extract_phenotype_vector(proteome_bed,phenotype_id)

if (length(transcriptomic_dat) > 0){
    
resid_transcriptomic <- transcriptomic_dat %>% residualize_molecular_phenotype_data(transcriptome_covar)
resid_proteomic <- proteomic_dat %>% residualize_molecular_phenotype_data(proteome_covar)

    
    
sample_ids <- proteome_bed %>% select(-1,-2,-3,-4) %>% colnames()
genotype_data <- proteome_bed %>%  
        extract_genotype_vector(phenotype_id,VCF_path) %>% 
        select(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,all_of(sample_ids))

genotype_matrix <- genotype_data %>% clean_genotype_data
phenotype_data <- list(outcome_1 =resid_transcriptomic,outcome_2 = resid_proteomic )

try(res <- colocboost(X = genotype_matrix[ , colSums(is.na(genotype_matrix)) == 0],
           Y = phenotype_data
           )
    )
try(res$name <- phenotype_id)
res
    }
}


