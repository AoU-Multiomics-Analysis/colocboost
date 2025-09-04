task colocboost {

    File GenotypeDosage 
    File GenotypeDosageIndex 
    File BedFile 
    File Covars
    Array[File] SumstatsGWAS
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
        --SumstatsGWAS ${SumstatsGWAS} \
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
        File colocboost_res = "~{phenotype_id}_colocboost_res.RDS"
    }
}


workflow colocboost_wdl {
    call colocboost
}
