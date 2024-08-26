version 1.0

task FilterVCF {
    input {
        File vcf_file
        String chrom
        String basename
        File indexed_vcf  # Add this input to take in the .tbi file
    }

    Int disk_size = 1 + 5*ceil(size([vcf_file, indexed_vcf], "GB"))

    command {
        echo 'FilterVCF'
        echo $vcf_file
        echo ${chrom}
        echo ${basename}
        echo ${indexed_vcf}
        bcftools view -r ${chrom} ${vcf_file} | \
        bcftools filter -e 'TYPE="indel"' -Oz -o ${basename}.${chrom}.vcf.gz
    }

    output {
        File filtered_vcf = "${basename}.${chrom}.vcf.gz"
    }
    
    runtime {
        cpu:                    2
        memory:                 8 + " GiB"
        disks: "local-disk " +  disk_size + " SSD" 
        bootDiskSizeGb:        10
        preemptible:            0
        maxRetries:             0
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task IndexFilteredVCF {
    input {
        File filtered_vcf
    }

    Int disk_size = 1 + 5*ceil(size([filtered_vcf], "GB"))

    command {
        bcftools index -t ${filtered_vcf}
    }

    output {
        File indexed_filtered_vcf = "${filtered_vcf}.tbi"
    }

    runtime {
        cpu:                    2
        memory:                 8 + " GiB"
        disks: "local-disk " +  disk_size + " SSD" 
        bootDiskSizeGb:        10
        preemptible:            0
        maxRetries:             0
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ConvertToZarr {
    input {
        File filtered_vcf
        String basename
        String chrom
    }

    Int disk_size = 1 + 5*ceil(size([filtered_vcf], "GB"))
    
    
    command <<<
    echo 'initializing vcf_to_zarr'
    echo ${filtered_vcf}
    echo ${basename}.${chrom}.zarr
        python3 <<CODE
import allel
input_file_py = f'~{filtered_vcf}'
zarr_output = f'~{basename}.~{chrom}.zarr'
print('test', input_file_py, zarr_output)
allel.vcf_to_zarr(input_file_py, 
                  zarr_output, 
                  fields='*', 
                  overwrite=False)
CODE

        tar czf ${basename}.${chrom}.zarr.tar.gz ${basename}.${chrom}.zarr
    >>>
    output {
        File zarr_tar = "${basename}.${chrom}.zarr.tar.gz"
    }

    runtime {
        cpu:                    2
        memory:                 8 + " GiB"
        disks: "local-disk " +  disk_size + " SSD" 
        bootDiskSizeGb:        10
        preemptible:            0
        maxRetries:             0
        docker:                 "mcfonsecalab/variantutils:0.8"
    }
}

workflow ConvertVCFtoZarr {
    input {
        File vcf_file
        String basename
        File indexed_vcf  # Add this input to take in the .tbi file
        String chrom      # Add this input for the chromosome
    }

    call FilterVCF { input: vcf_file = vcf_file, chrom = chrom, basename = basename, indexed_vcf = indexed_vcf }
    call IndexFilteredVCF { input: filtered_vcf = FilterVCF.filtered_vcf }
    call ConvertToZarr { input: filtered_vcf = FilterVCF.filtered_vcf, basename = basename, chrom = chrom }

    output {
        File zarr_tar = ConvertToZarr.zarr_tar
    }
}
