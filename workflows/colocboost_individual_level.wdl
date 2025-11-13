version 1.0 


task colocboost {
    input {
        File PhenotypeTable
        File ExpressionMatrix
        File SplicingMatrix
        File ProteinMatrix
        File ExpressionCovars
        File SplicingCovars
        File ProteinCovars
        File PlinkBedGenotypes
        File SampleList
        String OutputPrefix
    }
    
    command <<<
    Rscript /tmp/colocboost_v2.R \
        --PhenotypeTable ~{PhenotypeTable} \
        --ExpressionMatrix ~{ExpressionMatrix} \
        --SplicingMatrix ~{SplicingMatrix} \
        --ProteinMatrix ~{ProteinMatrix} \
        --ExpressionCovars ~{ExpressionCovars} \
        --SplicingCovars ~{SplicingCovars} \
        --ProteinCovars ~{ProteinCovars} \
        --PlinkBedGenotypes ~{PlinkBedGenotypes} \
        --SampleList ~{SampleList} \
        --OutputPrefix ~{OutputPrefix}

    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/colocboost:main"
        cpu: 1 
        memory: " 16 GB" 
        disks: "local-disk 200 HDD"
    }
    
    output {
        File ColocRes = "~{OutputPrefix}_colocboost_res.rds"
    }
}



workflow ColocboostWorkflow {
    input {
        File PhenotypeTable
        File ExpressionMatrix
        File SplicingMatrix
        File ProteinMatrix
        File ExpressionCovars
        File SplicingCovars
        File ProteinCovars
        File PlinkBedGenotypes
        File SampleList
        String OutputPrefix
    }
    
    call colocboost {
        input:
            PhenotypeTable = PhenotypeTable,
            ExpressionMatrix = ExpressionMatrix,
            SplicingMatrix = SplicingMatrix,
            ProteinMatrix = ProteinMatrix,
            ExpressionCovars = ExpressionCovars,
            SplicingCovars = SplicingCovars,
            ProteinCovars = ProteinCovars,
            PlinkBedGenotypes = PlinkBedGenotypes,
            SampleList = SampleList,
            OutputPrefix = OutputPrefix
    }
    output {
        File ColocSummary = colocboost.ColocRes 
    }
}
