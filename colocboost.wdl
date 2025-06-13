version 1.0

task stream_vcf {
    input {
        File VCF
        File VCF_index
        Int disk_space
    }
    
    parameter_meta {
        VCF: {
          description: "Cloud VCF file",
          localization_optional: true
        }
        VCF_index: {
          description: "Cloud VCF index",
          localization_optional: true
        }
    }

    command {
        bcftools view ${VCF} -r chr1:1000171-1000172
    }

    runtime {
        # Specify the resources required for the task
        docker: "quay.io/biocontainers/bcftools:1.22--h3a4d415_0"  # Example Docker image for tabix
        cpu: 1  # Default CPU allocation
        memory: "4 GB"  # Default memory allocation
        disks: "local-disk ${disk_space} HDD"
    }

    output {
        Array[File] out_vcfs = glob("*.vcf.gz")  # Collect all output VCF files
    }
}

task split_vcf {
    input {
        File VCF
        File VCF_index
        File proteome_bed
        Int padding
        Int disk_space
    }

    command {
        echo ${proteome_bed}

        #bgzip -d --output proteome.bed "${proteome_bed}"

        while IFS=$'\t' read -r chr start_pos end_pos name; do
            new_start=$((start_pos - padding))
            new_end=$((end_pos + padding))
            if (( new_start < 1 )); then new_start=1; fi

            region="$chr:$new_start-$new_end"
            out_vcf="$name.vcf.gz"

            echo $new_start
            echo $new_end
            echo $region
            echo $out_vcf

            # Extract the region from the VCF using tabix
            bcftools view -r $region -Oz "${VCF}"  > "$out_vcf"
            bcftools index -t "$out_vcf"  # Index the output VCF
        done < ${proteome_bed}
    }

    runtime {
        # Specify the resources required for the task
        #docker: "quay.io/biocontainers/tabix:0.2.5--0"  # Example Docker image for tabix
        docker: "epadhi/bcftools_bgzip:latest"
        cpu: 1  # Default CPU allocation
        memory: "4 GB"  # Default memory allocation
        disks: "local-disk ${disk_space} HDD"
    }

    output {
        Array[File] out_vcfs = glob("*.vcf.gz")  # Collect all output VCF files
        Array[File] out_index = glob("*.vcf.gz.tbi")  # Collect all output VCF files
    }
}



task colocboost {

    #File VCF
    #File VCF_index
    input {
        File transcriptome_bed 
        File proteome_bed
        File transcriptome_covars 
        File proteome_covars
        String phenotype_id
    
        String docker_image
        Int memory 
        Int disk_space 
        Int num_threads
    }

    command{
    echo "/src"
    ls /src
    echo "workdir"
    ls .
    cp /src/* .
    Rscript run_colocboost.R \
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
    call stream_vcf
    call split_vcf
    call colocboost
}
