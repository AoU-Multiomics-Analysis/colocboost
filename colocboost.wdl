task colocboost {

    #File VCF
    #File VCF_index
    File transcriptome_bed 
    File proteome_bed
    File transcriptome_covars 
    File proteome_covars
    String phenotype_id

    String docker_image
    Int memory 
    Int disk_space 
    Int num_threads

    command{
    ls
    ls /
    ls /cromwell_root
    tree .

    VCF_URL=$(gcloud storage sign-url gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/QTL_5500_subset/VCF/aou_rnaseq_5500.filtered.vcf.bgz --duration=1h --output-url)
    TBI_URL=$(gcloud storage sign-url gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/QTL_5500_subset/VCF/aou_rnaseq_5500.filtered.vcf.bgz.tbi --duration=1h --output-url)
    bcftools view "$VCF_URL" -r chr1:1000171-chr1:1000172

    Rscript run_colocboost.R \
    #    --vcf ${VCF} \
        --transcriptome_bed ${transcriptome_bed} \
        --proteome_bed ${proteome_bed}  \
        --transcriptome_covars ${transcriptome_covars}  \
        --proteome_covars ${proteome_covars} \
        --phenotype_id ${phenotype_id} 
    }
    runtime {
    docker: docker_image
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    bootDiskSizeGb: 25
    cpu: num_threads
    }
    
    output {
    
    # i dont think this is right but the phenotype id should go in the 
    # output
    File colocboost_res = "~{phenotype_id}_colocboost_res.RDS"

    }
}


workflow colocboost_wdl {
    call colocboost
}
