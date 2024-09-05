version 1.0

task FilterVCF {
    input {
        File vcf_file
        String chrom
        String basename
        File indexed_vcf
    }

    Int disk_size = 1 + 5*ceil(size([vcf_file, indexed_vcf], "GB"))

    command {
        echo 'FilterVCF'
        echo $vcf_file
        echo ${chrom}
        echo ${basename}
        bcftools view -r ${chrom} ${vcf_file} -Oz -o ${basename}.${chrom}.vcf.gz
    }

    output {
        File filtered_vcf = "${basename}.${chrom}.vcf.gz"
    }
    
    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " SSD" 
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " SSD" 
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ConvertToZarr {
    input {
        File filtered_vcf
        String chrom
        String basename
    }

    Int disk_size = 1 + 5*ceil(size(filtered_vcf, "GB"))

    command <<<
        set -x
        # set -euxo pipefail
        python3 <<EOF
        import allel
        import zarr
        vcfs = "~{filtered_vcf}"
        target = "~{basename}.~{chrom}.zarr"
        allel.vcf_to_zarr(vcfs, target, fields = ['variants/CHROM', 'variants/POS', 'variants/numalt', 'calldata/AF'], alt_number = 100)
        EOF
        echo "Done converting to zarr."
        echo "Tarring output..."
        find .
        tar -cf ~{basename}.~{chrom}.zarr.tar ~{basename}.~{chrom}.zarr
    >>>

    output {
        File zarr_file = "~{basename}.~{chrom}.zarr.tar"
    }

    runtime {
        cpu: 2
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " SSD" 
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/allelimage:latest3"
    }
}

workflow ConvertVCFtoZarr {
    input {
        File vcf_file
        String basename
        File indexed_vcf
        String chrom
       # File filtered_vcf
    }

    call FilterVCF { 
        input: vcf_file = vcf_file, chrom = chrom, basename = basename, indexed_vcf = indexed_vcf 
        }
   # call IndexFilteredVCF { 
       #input: filtered_vcf = FilterVCF.filtered_vcf }
    call ConvertToZarr { input: filtered_vcf = FilterVCF.filtered_vcf, basename = basename, chrom = chrom }
   # call ConvertToZarr { input: filtered_vcf = filtered_vcf, basename = basename, chrom = chrom }
  # call ConvertToZarr { input: filtered_vcf = filtered_vcf }

    output {
        File zarr_file = ConvertToZarr.zarr_file
    }
}