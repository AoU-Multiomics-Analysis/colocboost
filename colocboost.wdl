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
    echo "/src"
    ls /src
    echo "workdir"
    ls .
    Rscript /src/run_colocboost.R \
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
