version 1.0

task split_vcf_by_chromosome {
    input {
        File VCF
        File VCF_index
        Int disk_space
    }

    command <<<
        # Check if the input VCF file is valid
        bcftools view -h ${VCF} > /dev/null

        # Extract the list of chromosomes from the VCF file
        bcftools view -h ${VCF} | grep "^##contig" | sed 's/##contig=<ID=//;s/,.*//' > chromosomes.txt

        # Split the VCF file by chromosome
        while read chr; do
            out_vcf="${chr}.vcf.bgz"
            bcftools view -r $chr -Oz ${VCF} > $out_vcf
            bcftools index -t $out_vcf
        done < chromosomes.txt
    >>>

    runtime {
        docker: "epadhi/bcftools_bgzip:latest"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ${disk_space} HDD"
    }

    output {
        Array[File] chrom_vcfs = glob("*.vcf.gz")
        Array[File] chrom_indexes = glob("*.vcf.gz.tbi")
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

    command <<<
        echo ${proteome_bed}
        echo "Padding"
        echo ${padding}

        while IFS=$'\t' read -r chr start_pos end_pos name; do
            new_start=$((start_pos - ${padding}))
            new_end=$((end_pos + ${padding}))
            if (( new_start < 1 )); then new_start=1; fi

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
            bcftools view -r $region -Oz "${VCF}" > "$out_vcf"
            bcftools index -t "$out_vcf"  # Index the output VCF
        done < ${proteome_bed}
    >>>

    runtime {
        docker: "epadhi/bcftools_bgzip:latest"
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
        Int disk_space 
        Int num_threads
    }

    command <<<
        cp /src/* .
        Rscript run_colocboost.R \
            --transcriptome_bed ${transcriptome_bed} \
            --proteome_bed ${proteome_bed}  \
            --transcriptome_covars ${transcriptome_covars}  \
            --proteome_covars ${proteome_covars} \
            --phenotype_id ${phenotype_id} 
    >>>

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
    input {
        File VCF
        File VCF_index
        File transcriptome_bed
        File proteome_bed
        File transcriptome_covars
        File proteome_covars

        String docker_image
        Int memory
        Int disk_space
        Int num_threads
    }

    call split_vcf_by_chromosome {
        input:
            VCF = VCF,
            VCF_index = VCF_index,
            disk_space = disk_space
    }

    scatter (chrom_index in range(length(split_vcf_by_chromosome.chrom_vcfs))) {
        call split_vcf {
            input:
                VCF = split_vcf_by_chromosome.chrom_vcfs[chrom_index],
                VCF_index = split_vcf_by_chromosome.chrom_indexes[chrom_index],
                proteome_bed = proteome_bed,
                padding = 1000,  # Adjust padding as needed
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
