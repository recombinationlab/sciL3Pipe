
def log_parse(log):
    '''
    Parse split_hybrid.py log file(s).
    '''
    linecount = 0
    name, bam_1, bam_2  = [None] * 3
    for line in log:
        linecount += 1

        if line.rstrip().endswith('log'):
            name = line.rstrip()
        elif line.startswith('BAM 1'):
            bam_1 = line.rstrip().split(':')[-1]
        elif line.startswith('BAM 2'):
            bam_2 = line.rstrip().split(':')[-1]
            yield name, bam_1, bam_2
            name, bam_1, bam_2  = [None] * 3
        else:
            continue


if config["run_split"]:
    rule split_hybrid:
        '''
        Assume the BAM files are coordinate sorted

        Split is done based on ensembl (1,2,...,X,Y - human) vs ucsc (chr1,chr2,...,chrX,chrY - mouse) coordinates
        Is based on how the hybrid reference was made
        '''
        input:
            # lambda wildcards: FILES[os.path.basename(wildcards.sample)]
            out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam",
            out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam.bai"
        output:
            out_dir + f"split_bam_{aligner}_{assembly}_hybrid/{{sample}}/{{sss}}/{prefix_1}/{{i}}.{prefix_1}.bam",
            out_dir + f"split_bam_{aligner}_{assembly}_hybrid/{{sample}}/{{sss}}/{prefix_2}/{{i}}.{prefix_2}.bam",
            out_dir + f"logs/{{sample}}/{{sss}}/{{i}}_hybrid_split.log"
        conda:
            "../envs/sciL3_env.yaml"
        log:
            out_dir + f"logs/{{sample}}/{{sss}}/{{i}}_hybrid_split.err"
        params:
            output_dir = out_dir + f"split_bam_{aligner}_{assembly}_hybrid//{{sample}}/{{sss}}/"
        group:
            "groupjob_hyb_split"
        shell:
            '''
            python workflow/scripts/split_hybrid.py \
            -i {input[0]} \
            -o {params.output_dir} \
            -1 {prefix_1} \
            -2 {prefix_2} \
            {bam_filter} > {output[2]} 2> {log}
            '''


    def input_qc(wildcards):
        if sss_bam_start:
            log_files = []
            for s in SSS_FILES.keys():
                if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                    sss_keys = list(SSS_FILES.get(s).keys())
                    for i in sss_keys:
                        # print('i', i)
                        # log_files.extend([out_dir + f"merged_sss_bam/{s}/{sss}.PE.{aligner}.{assembly}.srt.bam"])
                        output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                        # print('output_j', output_j)
                        outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                        # print('outputs_j', outputs_j)
                        for j in outputs_j:
                            log_files.extend(expand(out_dir + f"logs/{s}/{i}/{j}_hybrid_split.log"))
            # print("summary:", log_files)
        elif trimmed_start:
            log_files = []
            for s in TRM_FILES.keys():
                if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                    sss_keys = list(TRM_FILES.get(s).keys())
                    for i in sss_keys:
                        output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                        outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                        for j in outputs_j:
                            log_files.extend(expand(out_dir + f"logs/{s}/{i}/{j}_hybrid_split.log"))

        else:

            outputs_i_wc = glob.glob(f"{checkpoints.split_pe.get(**wildcards).output}/*R1*.fastq.gz")
            outputs_i_wc = [output.split('/')[-1].split('.')[2] for output in outputs_i_wc]
            # print('outputs_i_wc', outputs_i_wc)

            log_files = []
            for i_wc in outputs_i_wc:
                # print('i_wc', i_wc)
                outputs_i = glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards, sfq=i_wc).output}/*R1.ordered.*")
                outputs_i = [output.split('/')[-1].split('.')[0] for output in outputs_i]
                
                # print('outputs_i', outputs_i)
                for i in outputs_i:
                    output_j = glob.glob(f"{checkpoints.split_bam.get(**wildcards, sss=i).output}/*.bam")
                    outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]

                    for j in outputs_j:
                        log_files.extend(expand(out_dir + f"logs/{{sample}}/{i}/{j}_hybrid_split.log", sample=wildcards.sample))
                        
            # print("summary:", log_files)

        return log_files


    rule qc:
        '''
        Only output subset of BAM files, based on human mouse mix
        '''
        input:
            # [OUTPUT_DIR[s] + f"logs/{s}_hybrid_split.log" for s in ALL_SAMPLES]
            # out_dir + f"logs/{{sample}}/{{sss}}/{{i}}_hybrid_split.log"
            input_qc
        output:
            out_dir + "logs/{sample}_qc.tsv"
        run:
            rows = []
            for f in input:
                with open(f) as tsv:
                    for name, bam_1, bam_2 in log_parse(tsv):
                        total = int(bam_1) + int(bam_2)
                        bam_1_pot = (int(bam_1)/total)*100
                        bam_2_pot = (int(bam_2)/total)*100
                        rows.append([f, bam_1, bam_2, bam_1_pot, bam_2_pot, total])

            df = pd.DataFrame(rows, columns=['Files', 'BAM_1', 'BAM_2', 'BAM_1_PofT', 'BAM_2_PofT', 'Total'])
            df.to_csv(output[0], sep='\t', index=False)



rule filter_bam:
    '''
    Filter excessive soft-clipping
    '''
    input:
        # f"{{output_dir}}{process_only}/{{sample}}.{process_only}.bam" if config["run_split"] else lambda wildcards: FILES[os.path.basename(wildcards.sample)]
        out_dir + f"split_bam_{aligner}_{assembly}_hybrid/{{sample}}/{{sss}}/{process_only}/{{i}}.{process_only}.bam" if config["run_split"] else out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam"
    output:
        # f"{{output_dir}}{process_only}/filtered/{{sample}}.{process_only}.filt.bam"
        out_dir + f"split_bam_{aligner}_{assembly}_filtered/{{sample}}/{{sss}}/{process_only}/{{i}}.{process_only}.filt.bam"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + "logs/{sample}/{sss}.{i}_filtered.log"
    threads:
        1
    group:
        "groupjob_filt_bam"
    shell:
        '''
        python workflow/scripts/filter_bam.py \
        -i {input} \
        -o {output} \
        --threads {threads} \
       --paired \
       --insert_max 2000 \
        --sort \
        --no_plot \
        --max_ratio 0.5 \
        --filter_alt \
        {bam_filter_2} &> {log}
        '''




def filtered_summary_input(wildcards):
    if sss_bam_start:
        filt_files = []
        for s in SSS_FILES.keys():
            if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                sss_keys = list(SSS_FILES.get(s).keys())
                for i in sss_keys:
                    # print('i', i)
                    # filt_files.extend([out_dir + f"merged_sss_bam/{s}/{sss}.PE.{aligner}.{assembly}.srt.bam"])
                    output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                    # print('output_j', output_j)
                    outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                    # print('outputs_j', outputs_j)
                    for j in outputs_j:
                        filt_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}_filtered/{s}/{i}/{process_only}/{j}.{process_only}.filt.bam"))
        # print("summary:", filt_files)
    elif trimmed_start:
        filt_files = []
        for s in TRM_FILES.keys():
            if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                sss_keys = list(TRM_FILES.get(s).keys())
                for i in sss_keys:
                    output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                    outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                    for j in outputs_j:
                        filt_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}_filtered/{s}/{i}/{process_only}/{j}.{process_only}.filt.bam"))

    else:

        outputs_i_wc = glob.glob(f"{checkpoints.split_pe.get(**wildcards).output}/*R1*.fastq.gz")
        outputs_i_wc = [output.split('/')[-1].split('.')[2] for output in outputs_i_wc]
        # print('outputs_i_wc', outputs_i_wc)

        filt_files = []
        for i_wc in outputs_i_wc:
            # print('i_wc', i_wc)
            outputs_i = glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards, sfq=i_wc).output}/*R1.ordered.*")
            outputs_i = [output.split('/')[-1].split('.')[0] for output in outputs_i]
            
            # print('outputs_i', outputs_i)
            for i in outputs_i:
                output_j = glob.glob(f"{checkpoints.split_bam.get(**wildcards, sss=i).output}/*.bam")
                outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]

                for j in outputs_j:
                    filt_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}_filtered/{{sample}}/{i}/{process_only}/{j}.{process_only}.filt.bam",
                                              sample=wildcards.sample))

        # print("summary:", filt_files)

    return filt_files


rule filtered_summary_list:
    '''
    If split_hybrid, then BAM files that have fewer than 30% of total reads or less than 10000 reads will be excluded
    '''
    input:
        fsi=filtered_summary_input,
        qc=out_dir + "logs/{sample}_qc.tsv" if config["run_split"] else []
    output:
        out_dir + f"files_{aligner}_{assembly}/{{sample}}_{process_only}.filt.list"
    run:
        if config["run_split"]:
            qc = pd.read_csv(input.qc, sep="\t")
            
            # BAM has fewer than 30% of total reads or less than 10000 reads and will be excluded
            if process_only == 'mouse': # BAM_1
                samples = qc[(pd.to_numeric(qc['BAM_2_PofT']) < 30) & (pd.to_numeric(qc['Total']) > pd.to_numeric(config['cutoff']))]['Files']
            elif process_only == 'human': # BAM_2
                samples = qc[(pd.to_numeric(qc['BAM_1_PofT']) < 30) & (pd.to_numeric(qc['Total']) > pd.to_numeric(config['cutoff']))]['Files']

            samples_pass = []
            for s in samples:
                samples_pass.append(os.path.basename(s).split("_hybrid_split.log")[0])
          
            filt_input = [i for i in input.fsi if os.path.basename(i).split(f".{process_only}.filt.bam")[0] in samples_pass]

            with open(output[0], 'w') as out:
                out.write('\n'.join(filt_input))
            
        else:
            with open(output[0], 'w') as out:
                out.write('\n'.join(input.fsi))

# rule idxstats_qc:
#     '''
#     Summary of idxstats for all BAM files
#     '''
#     input:
#         [OUTPUT_DIR[s] + f"{process_only}/{s}.{process_only}.bam" for s in ALL_SAMPLES] if config["run_split"] else list(chain(*[FILES[s] for s in ALL_SAMPLES]))
#     output:
#         qc_folder + "logs/idxstats_summary.tsv"
#     conda:
#         "../envs/sciL3_env.yaml"
#     params:
#         dir = ' '.join(config["input_folder"]),
#         format = config["style"]
#         # format = 'ensembl' if process_only == 'human' else 'UCSC'
#     shell:
#         '''
#         echo {params.dir}
#         python scripts/aggregate_idxstats.py \
#         -d {params.dir} \
#         -o {output} \
#         --format {params.format} \
#         --regex "*.bam"
#         '''
