task colocboost {

    File GenotypeDosage 
    File GenotypeDosageIndex 
    File BedFile 
    File Covars 
    String PhenotypeID 

    String docker_image
    Int memory 
    Int disk_space 
    Int num_threads

    command{
    Rscript colocboost_summarystats.R \
        --GenotypeDosage ${GenotypeDosage} \
        --BedFile ${BedFile} \
        --Covars ${Covars}  \
        --PhenotypeID ${PhenotypeID} 
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
