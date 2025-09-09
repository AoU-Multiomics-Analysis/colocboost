
task colocboost {

    File GenotypeDosage 
    File GenotypeDosageIndex 
    File BedFile 
    File Covars
    Array[File] SumstatsGWAS
    String PhenotypeID 
    
    Int NumPrempt
    Int memory 
    Int disk_space 
    Int num_threads

    command{
    Rscript /tmp/colocboost_summarystats.R \
        --GenotypeDosage ${GenotypeDosage} \
        --BedFile ${BedFile} \
        --Covars ${Covars}  \
        --SumstatsGWAS ${sep=" " SumstatsGWAS} \
        --PhenotypeID ${PhenotypeID} 
    }
    runtime {
    docker: "ghcr.io/aou-multiomics-analysis/colocboost:main"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    bootDiskSizeGb: 25
    preemptible: "${NumPrempt}"
    cpu: num_threads
    }
    
    output {
        File colocboost_res = "~{PhenotypeID}_colocboost_res.RDS"
    }
}


workflow colocboost_wdl {
    call colocboost
}
