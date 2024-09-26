version 1.0

task Createh5 {
    input {
        Array[File] standard3000_files
    }

    Int disk_size = 1 + 5 * ceil(size(standard3000_files, "GB"))

    command <<<

        python3 <<EOF
        import numpy as np
        import h5py

        # Array of file paths
        file_paths = ~{sep=' ' standard3000_files}
        combined_file = 'combined_14_files.h5'

        combined_samples = []  # To store samples from all files
        other_keys_data = {}  # To temporarily store other keys' data

        for file_path in file_paths.split():
            with h5py.File(file_path, 'r') as f:
                for key in f.keys():
                    if key == 'samples':
                        # Concatenate 'samples' datasets
                        samples = f['samples'][:]
                        combined_samples.append(samples)
                    else:
                        # Store the first occurrence of each key's data
                        if key not in other_keys_data:
                            other_keys_data[key] = f[key][:]

        # Combine all 'samples' from the files
        combined_samples = np.concatenate(combined_samples, axis=0)

        with h5py.File(combined_file, 'w') as combined_h5:
            combined_h5.create_dataset('samples', data=combined_samples)
            for key, data in other_keys_data.items():
                combined_h5.create_dataset(key, data=data)

        EOF
    >>>

    output {
        File standard3000 = "combined_14_files.h5"
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/1000-3000-comparison:latest"
    }
}

workflow CompileData {
    input {
        Array[File] standard3000_files
    }

    call Createh5 { input: 
            standard3000_files = standard3000_files
    }

    output {
        File standard3000 = Createh5.standard3000
    }
}
