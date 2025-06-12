task colocboost {

    File VCF
    File transcriptome_bed 
    File proteome_bed
    File transcriptome_covars 
    File proteome_covars 
    String phenotype_id
    
    command{
    Rscript run_colocboost.R \
        --vcf \
        --transcriptome_bed \
        --proteome_bed \
        --transcriptome_covars \
        --proteome_covars \
        --phenotype_id \
        --output 
    }
    runtime {
    docker: ""
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    bootDiskSizeGb: 25
    cpu: "${num_threads}"
    preemptible: "${num_preempt}"
    }
    
    output {
    
    # i dont think this is right but the phenotype id should go in the 
    # output
    File colocboost_res = "${phenotype_id}_colocboost_res.RDS"

    }
}


workflow colocboost {
    call colocboost
}
