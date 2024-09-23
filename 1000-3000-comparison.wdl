version 1.0

task FilterFile {
    input {
        File input_file
        File pos_mask_dict
        File sample_mask_dict
        String interest
    }

    Int disk_size = 1 + 5*ceil(size([input_file, pos_mask_dict, sample_mask_dict], "GB"))

    command <<<

        python3 <<EOF
        import numpy as np 
        import h5py
        import json

        read_file = h5py.File("~{input_file}", 'r')

        with open("~{pos_mask_dict}", 'r') as f:
            pos_mask_dict = json.load(f)

        with open("~{sample_mask_dict}", 'r') as f:
            sample_mask_dict = json.load(f)

        final_data_sets = {}

        def make_datasets(dataset, comparison_of_interest):
            subset_original = {"genotypes": [], "ad": [], "ab": [], "position": []}
            start_idx = 0
            for chromosome in range(1, 15):
                array_length = len(dataset['POS_' + str(chromosome)][::])
                end_idx = start_idx + array_length
                pos_mask = pos_mask_dict[comparison_of_interest]
                sample_mask = np.asarray(sample_mask_dict[comparison_of_interest])
                subset_pos_mask = np.asarray(pos_mask[start_idx:end_idx])
                start_idx = end_idx
                pos_array = np.asarray(dataset["POS_" + str(chromosome)][::])
                avail_pos = pos_array[~subset_pos_mask]
                
                for pos in avail_pos:
                    subset_original['position'].append(str(chromosome) + ":" + str(pos))
                
                genotypes = np.asarray(dataset[str(chromosome)][::].T)
                ab = np.asarray(dataset['p_intra_' + str(chromosome)][::])
                ad = np.asarray(dataset["DP_" + str(chromosome)][::])
                ad_pos_filter = ad[~sample_mask]
                subset_original["ad"].append(ad_pos_filter.T[~subset_pos_mask])
                ab_pos_filter = ab[~sample_mask]
                subset_original["ab"].append(ab_pos_filter.T[~subset_pos_mask])
                gen_pos_filter = genotypes[~sample_mask]
                subset_original["genotypes"].append(gen_pos_filter.T[~subset_pos_mask])
            
            subset_final = {"genotypes": [], "ad": [], "ab": [], 'position': []}
            original_samples = np.asarray([x.decode('UTF-8') for x in dataset["samples"][::]])
            subset_final['samples'] = original_samples[~sample_mask]
            subset_final['genotypes'] = np.concatenate(subset_original['genotypes'])
            subset_final['ad'] = np.concatenate(subset_original['ad'])
            subset_final['ab'] = np.concatenate(subset_original['ab'])
            subset_final['position'] = subset_original['position']
            final_data_sets[comparison_of_interest] = subset_final
        
        make_datasets(read_file, "~{interest}")

        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)
                
        with open("~{interest}.json", "w") as outfile:
            json.dump(final_data_sets["~{interest}"], outfile, cls=NumpyEncoder)

        EOF
    >>>

    output {
        File new_file = "${interest}.json"
    }
    
    runtime {
        cpu: 2
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/1000-3000-comparison:latest"
    }
}

workflow SubsetData {
    input {
       File input_file
       File pos_mask_dict
       File sample_mask_dict
       String interest
    }

   call FilterFile { input: input_file = input_file, pos_mask_dict = pos_mask_dict, sample_mask_dict = sample_mask_dict, interest = "interest" } # You might want to provide interest dynamically

    output {
        File new_file = FilterFile.new_file
    }
}
