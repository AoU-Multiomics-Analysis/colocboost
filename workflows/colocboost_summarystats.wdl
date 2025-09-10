version 1.0
workflow colocboost_wdl {
    input { 
        File GenotypeDosage 
        File GenotypeDosageIndex 
        File BedFile 
        File Covars
        Array[File] SumstatsGWAS
        Array[File] SumstatsGWASIndex
        String PhenotypeID 
        Boolean Partialize        
        Int NumPrempt
        Int memory 
        Int disk_space 
        Int num_threads

    }

    call colocboost {
        input:
            GenotypeDosage = GenotypeDosage,
            GenotypeDosageIndex = GenotypeDosageIndex,
            BedFile = BedFile,
            Covars = Covars,
            SumstatsGWAS = SumstatsGWAS,
            SumstatsGWASIndex = SumstatsGWASIndex,
            PhenotypeID = PhenotypeID,
            NumPrempt = NumPrempt,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            Partialize  = Partialize 

    }
    output {
        File colocboost_res = colocboost.colocboost_res
        File colocboost_summary = colocboost.colocboost_summary
    }
}

    task colocboost {
        input {
            File GenotypeDosage 
            File GenotypeDosageIndex 
            File BedFile 
            File Covars
            Array[File] SumstatsGWAS
            Array[File] SumstatsGWASIndex
            String PhenotypeID 
            Boolean Partialize        
      
            Int NumPrempt
            Int memory 
            Int disk_space 
            Int num_threads
        }
        command <<<
        Rscript /tmp/colocboost_summarystats.R \
            --GenotypeDosage ~{GenotypeDosage} \
            --BedFile ~{BedFile} \
            --Covars ~{Covars}  \
            --SumstatsGWAS ~{sep=" " SumstatsGWAS} \
            --PhenotypeID ~{PhenotypeID} \
            --Partialize ~{Partialize}
        >>>
        runtime {
            docker: "ghcr.io/aou-multiomics-analysis/colocboost:main"
            memory: "~{memory}GB"
            disks: "local-disk ~{disk_space} HDD"
            bootDiskSizeGb: 25
            preemptible: "~{NumPrempt}"
            cpu: num_threads
        }
        
        output {
            File colocboost_res = "~{PhenotypeID}_colocboost_res.RDS"
            File colocboost_summary = "~{PhenotypeID}_colocboost_summary.RDS"

        }
    }


