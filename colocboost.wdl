task colocboost {

    File VCF
    File transcriptome_bed 
    File proteome_bed
    File transcriptome_covars 
    File proteome_covars 
    String phenotype_id
    
    command{



    }
    runtime {


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
