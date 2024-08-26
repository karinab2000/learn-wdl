version 1.0

task ConvertToZarr {
    input {
        File filtered_vcf
        String basename
        String chrom
    }

    Int disk_size = 1 + 5*ceil(size([filtered_vcf], "GB"))

    command <<<
        echo 'initializing vcf_to_zarr'
        echo '~{filtered_vcf}'
        echo "~{basename}.~{chrom}.zarr"
        
        python3 <<CODE
        import allel
        #import os


        input_file_py = "~{filtered_vcf}"
        #input_dir_path = os.path.dirname(os.path.realpath(__file__))
        zarr_output = "~{basename}.~{chrom}.zarr"
        print("test", input_file_py, zarr_output)
        allel.vcf_to_zarr("~{filtered_vcf}", 
                        "~{basename}.~{chrom}.zarr", 
                        fields="*", 
                        overwrite=False)
        CODE

        #tar czf ${basename}.${chrom}.zarr.tar.gz ${basename}.${chrom}.zarr
    >>>

    output {
        String zarr_tar = 'test'
        #File zarr_tar = "~{basename}.~{chrom}.zarr.tar.gz"
    }

    runtime {
        cpu:                    2
        memory:                 20 + " GiB"
        disks: "local-disk " +  disk_size + " SSD" 
        bootDiskSizeGb:        10
        preemptible:            0
        maxRetries:             0
        docker:                 "karinab2000/allelimage:latest"
    }
}


workflow ConvertVCFtoZarr {
    input {
        File vcf_file
        String basename
        File indexed_vcf  # Add this input to take in the .tbi file
        String chrom      # Add this input for the chromosome
        File filtered_vcf

    }

    #call FilterVCF { input: vcf_file = vcf_file, chrom = chrom, basename = basename, indexed_vcf = indexed_vcf }
    #call IndexFilteredVCF { input: filtered_vcf = FilterVCF.filtered_vcf }
    call ConvertToZarr { input: filtered_vcf = filtered_vcf, basename = basename, chrom = chrom }

    output {
        File zarr_tar = ConvertToZarr.zarr_tar
    }
}
