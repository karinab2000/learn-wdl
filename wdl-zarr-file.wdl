version 1.0

task FilterVCF {
    input {
        File vcf_file
        String chrom
        String basename
        File indexed_vcf  # Add this input to take in the .tbi file
    }

    command {
        bcftools view -r ${chrom} ${vcf_file} | \
        bcftools filter -e 'TYPE="indel"' -Oz -o ${basename}.${chrom}.vcf.gz
    }

    output {
        File filtered_vcf = "${basename}.${chrom}.vcf.gz"
    }
}

task IndexFilteredVCF {
    input {
        File filtered_vcf
    }

    command {
        bcftools index -t ${filtered_vcf}
    }

    output {
        File indexed_filtered_vcf = "${filtered_vcf}.tbi"
    }
}

task ConvertToZarr {
    input {
        File filtered_vcf
        String basename
        String chrom
    }

    command <<<
        python3 <<CODE
import allel

allel.vcf_to_zarr('${filtered_vcf}', 
                  '${basename}.${chrom}.zarr', 
                  fields='*', 
                  overwrite=False)
CODE
    >>>
    
    output {
        File zarr_output = "${basename}.${chrom}.zarr"
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
}
