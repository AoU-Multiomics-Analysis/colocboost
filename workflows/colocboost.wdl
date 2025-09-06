version 1.0

task split_vcf_by_chromosome {
    input {
        File VCF
        File VCF_index
        String chromosome
        Int disk_space
    }

    command <<<
        # Check if the input VCF file is valid
        echo ${VCF}
        echo "${VCF}"
        echo ~{VCF}
        echo "~{VCF}"
        echo "VCF file: ~{VCF}"
        echo "Checking file"
        bcftools view -h "~{VCF}" > /dev/null

        # Split the VCF file by chromosome
        echo "Processing chromosome: ~{chromosome}"
        out_vcf="~{chromosome}.vcf.bgz"
        echo "Output VCF: $out_vcf"
        bcftools view -r ~{chromosome} -Oz "~{VCF}" > $out_vcf
        bcftools index -t $out_vcf
    >>>

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.22--h3a4d415_0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ${disk_space} HDD"
    }

    output {
        File chrom_vcf = "~{chromosome}.vcf.bgz"
        File chrom_index = "~{chromosome}.vcf.bgz.tbi"
    }
}


task split_vcf {
    input {
        File VCF
        File VCF_index
        String chromosome
        File trimmed_proteome_bed
        Int padding
        Int disk_space
    }

    command <<<
        echo ~{trimmed_proteome_bed}
        echo "Padding"
        echo ~{padding}
        bcftools index -t ~{VCF}
        while IFS=$'\t' read -r chr start_pos end_pos name; do
            new_start=$((start_pos - ~{padding}))
            new_end=$((end_pos + ~{padding}))
            if (( new_start < 1 )); then new_start=1; fi
            if [[ "~{chromosome}" != "$chr" ]]; then continue; fi
            region="$chr:$new_start-$new_end"
            out_vcf="$name.vcf.gz"

            echo "Entry"
            echo $start_pos
            echo $end_pos
            echo "Padding Entry"
            echo $new_start
            echo $new_end
            echo $region
            echo $out_vcf

            # Extract the region from the VCF using tabix
            bcftools view -r $region -Oz "~{VCF}" > "$out_vcf"
            bcftools index -t "$out_vcf"  # Index the output VCF
        done < ~{trimmed_proteome_bed}
    >>>

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.22--h3a4d415_0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ${disk_space} HDD"
    }

    output {
        Array[File] out_vcfs = glob("*.vcf.gz")
        Array[File] out_index = glob("*.vcf.gz.tbi")
    }
}

task colocboost {
    input {
        File VCF
        File VCF_index
        File transcriptome_bed 
        File proteome_bed
        File transcriptome_covars 
        File proteome_covars
        String phenotype_id
    
        String docker_image
        Int memory
        String disk_type = "HDD"
        Int disk_space 
        Int boot_disk_space = 25
        Int num_threads
    }

    command <<<
        cp /src/* .
        tabix -p vcf ~{VCF}
        Rscript run_colocboost.R \
            --vcf ~{VCF} \
            --transcriptome_bed ~{transcriptome_bed} \
            --proteome_bed ~{proteome_bed}  \
            --transcriptome_covars ~{transcriptome_covars}  \
            --proteome_covars ~{proteome_covars} \
            --phenotype_id ~{phenotype_id} 
    >>>

    runtime {
        docker: docker_image
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} ${disk_type}"
        bootDiskSizeGb: boot_disk_space
        cpu: num_threads
    }
    
    output {
        File colocboost_res = "~{phenotype_id}_colocboost_res.RDS"
    }
}

workflow colocboost_wdl {
    input {
        File VCF_workflow
        File VCF_workflow_index
        File transcriptome_bed
        File trimmed_proteome_bed
        File proteome_bed
        File transcriptome_covars
        File proteome_covars

        String docker_image
        Int memory
        Int disk_space
        Int num_threads
    }

    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    scatter (chromosome in chromosomes) {
        call split_vcf_by_chromosome {
            input:
                VCF = VCF_workflow,
                VCF_index = VCF_workflow_index,
                chromosome = chromosome,
                disk_space = disk_space
        }

        call split_vcf {
            input:
                VCF = split_vcf_by_chromosome.chrom_vcf,
                VCF_index = split_vcf_by_chromosome.chrom_index,
                chromosome = chromosome,
                trimmed_proteome_bed = trimmed_proteome_bed,
                padding = 1000000,  # Adjust padding as needed
                disk_space = disk_space
        }

        scatter (index in range(length(split_vcf.out_vcfs))) {
            String vcf_file_name = basename(split_vcf.out_vcfs[index])
            String phenotype_id = sub(vcf_file_name, ".vcf.gz$", "")

            call colocboost {
                input:
                    VCF = split_vcf.out_vcfs[index],
                    VCF_index = split_vcf.out_index[index],
                    transcriptome_bed = transcriptome_bed,
                    proteome_bed = proteome_bed,
                    transcriptome_covars = transcriptome_covars,
                    proteome_covars = proteome_covars,
                    phenotype_id = phenotype_id,
                    docker_image = docker_image,
                    memory = memory,
                    disk_space = disk_space,
                    num_threads = num_threads
            }
        }
    }
}
