

def assign_sample(bam_path, sample_lib_sss, sample_tn5):
    '''
    Identify barcodes in file name

    Args:
        bam_path (str): Path to BAM file.
        sample_lib_sss (list[str]): List of lib_sss barcodes belonging to sample.
        sample_tn5 (list[str]): List of tn5 barcodes belonging to sample.
    '''
    # bam_pattern = re.compile(r"(?P<lib>[^_/]+)_(?P<sss>[ACGT]{6})\.(?P<tn5>[ACGT]{8})(?P<lig>[ACGT]{7})\.bam$")
    # bam_pattern = re.compile(r"(?P<lib_SSS>[^_/]+_[ACGT]{6})\.(?P<tn5>[ACGT]{8})(?P<lig>[ACGT]{7})\.bam$")
    bam_pattern = re.compile(r"(?P<lib_SSS>[^_/]+_[ACGT]{6})\.(?P<tn5>[ACGT]{8})(?P<lig>[ACGT]{7})(?P<suffix>(?:\.[^.]+)*)\.bam$")
    m = bam_pattern.search(bam_path)
    if not m:
        return None  # not matching pattern
    
    lib_sss = m.group("lib_SSS")
    tn5 = m.group("tn5")
    lig = m.group("lig")

    #  Only check if the list is non-empty
    # lib_ok = (lib in sample_lib) if sample_lib else True
    lib_sss_ok = (lib_sss in sample_lib_sss) if sample_lib_sss else True
    tn5_ok = (tn5 in sample_tn5) if sample_tn5 else True

    if lib_sss_ok and tn5_ok:
        return bam_path
    else:
        return None


def create_symlinks(paths, outdir):
    '''
    Create symlinks for a set of file paths into the output directory.
    
    Args:
        paths (set[str] or list[str]): Paths to the source files.
        outdir (str or Path): Path to the folder where symlinks should be created.
    '''
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)  # create output folder if needed

    for src_path in paths:
        src = Path(src_path)
        if not src.exists():
            print(f"Warning: {src} does not exist, skipping.")
            continue

        link = outdir / src.name  # symlink name same as source filename
        try:
            if link.exists():
                link.unlink()  # remove existing link/file
            link.symlink_to(src.resolve())  # create symlink
            print(f"Created symlink: {link} -> {src}")
        except Exception as e:
            print(f"Failed to create symlink for {src}: {e}")



        


def assign_sample_input(wildcards):
    # sample wildcard here is library (e.g. 'yi401')
    
    if sss_bam_start:
        split_files = list()
        for s in SSS_FILES.keys():
            if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                sss_keys = list(SSS_FILES.get(s).keys())
                for i in sss_keys:
                    output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                    outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                    for j in outputs_j:
                        split_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}/{s}/{i}/{j}.bam"))

        # print("summary:", split_files)
    elif trimmed_start:
        split_files = []
        for s in TRM_FILES.keys():
            if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                sss_keys = list(TRM_FILES.get(s).keys())
                for i in sss_keys:
                    output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                    outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                    for j in outputs_j:
                        split_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}/{s}/{i}/{j}.bam"))
                      
    else:
        outputs_i_wc = glob.glob(f"{checkpoints.split_pe.get(**wildcards).output}/*R1*.fastq.gz")
        outputs_i_wc = [output.split('/')[-1].split('.')[2] for output in outputs_i_wc]
        # print('outputs_i_wc', outputs_i_wc)

        split_files = list()
        for i_wc in outputs_i_wc:
            # print('i_wc', i_wc)
            outputs_i = glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards, sfq=i_wc).output}/*R1.ordered.*")
            outputs_i = [output.split('/')[-1].split('.')[0] for output in outputs_i]
            
            print('outputs_i', outputs_i)
            for i in outputs_i:
                output_j = glob.glob(f"{checkpoints.split_bam.get(**wildcards, sss=i).output}/*.bam")
                outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]

                for j in outputs_j:
                    split_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{i}/{j}.bam",
                                                sample=wildcards.sample))

    # print(' '.join(split_files))
    return set(split_files)





rule assign_sample:
    '''
    Using a dictionary of single cell BAMs, create sample folders with symlinks to sample BAMs
    '''
    input:
        out_dir + f"files_{aligner}_{assembly}/{{sample}}_{process_only}.filt.list" if config['run_filter'] else assign_sample_input
    output:
        out_dir + f"samples_assigned_{aligner}_{assembly}/{{sample}}.list"
    log:
        out_dir + f"logs/{{sample}}.{aligner}_{assembly}.samples_assigned.log"
    run:
        # if filtering run, use .list files as input
        if config['run_filter']:
            sample_files = list()
            for i in input:
                with open(i) as in_list:
                    for f in in_list:
                        sample_files.append(f.strip())

        # Parse split_file based on config
        with open(output[0], 'a') as out:
            for i, entry in enumerate(config['sample_assignment']):
                sample_to_assign = list(entry.keys())[0] 
                sample_info = list(entry.values())[0]  # gives {"lib_SSS": [...], "Tn5": [...]}
                sample_lib_sss = sample_info['lib_SSS']
                sample_tn5 = sample_info['Tn5']
                sample_lib = [i.split('_')[0] for i in sample_lib_sss]

                pass_files = set()
                if config['run_filter']:
                    for f in sample_files:
                        assignment = assign_sample(f, sample_lib_sss, sample_tn5)
                        if assignment is not None:
                            pass_files.add(f)
                else:
                    for f in input:
                        assignment = assign_sample(f, sample_lib_sss, sample_tn5)
                        if assignment is not None:
                            pass_files.add(f)
                    
                if wildcards.sample in sample_lib:
                    out.write(sample_to_assign + '\n' + '\n'.join(pass_files) + '\n')

                create_symlinks(pass_files, f"{out_dir}{sample_to_assign}/")

        print("All symlinks created successfully.")
        
        