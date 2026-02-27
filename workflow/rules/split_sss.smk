


checkpoint split_pe:
    '''
    Split large fastq files into smaller chunks that can be trimmed in parallel 
    '''
    wildcard_constraints:
        sample='|'.join(samples)
    input:
        r1 = out_dir + "merged/{sample}.R1.fastq.gz",
        r2 = out_dir + "merged/{sample}.R2.fastq.gz"
    output:
        temp(directory(out_dir + "fastqs_split/{sample}/"))
    log:
        stdout = out_dir + "logs/{sample}_fq_split.out",
        stderr = out_dir + "logs/{sample}_fq_split.err"
    params:
        lines=config['split_lines']
    conda:
        "../envs/sciL3_env.yaml"
    shell:
        ''' 
        seqkit split2 -1 {input.r1} -2 {input.r2} -s {params.lines} -O {output} -e .gz > {log.stdout} 2> {log.stderr}
        '''



    # def split_pe_aggregate(wildcards):
    #     checkpoint_output = checkpoints.split_pe.get(**wildcards).output[0]
    #     print('wc sample', wildcards.sample)
    #     out = expand(out_dir + "fastqs_split/{sample}/{i}.fastq.gz",
    #             sample=wildcards.sample,
    #             i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fastq.gz")).i)
    #     print('Split PE output:', out)
    #     return out


    # rule sample_json:
    #     input:
    #         split_pe_aggregate
    #     output:
    #         # touch(out_dir + f"logs/{{sample}}_{aligner}_{assembly}.align_all_done.txt")
    #         out_dir + "{sample}_split_samples.json"
    #     params:
    #         fastq_dir = out_dir + "fastqs_split/{sample}/"
    #     shell:
    #         '''
    #         python fastq2json.py --fastq_dir {params.fastq_dir} --split -o {output}
    #         '''


    # rule test:
    #     input:
    #         out_dir + "{sample}_split_samples.json"
    #     output:
    #         out_dir + "{sample}_test.json"
    #     run:
    #         with open(input[0], "r") as f:
    #             # Load the JSON data
    #             data = json.load(f)
    #             print(data)
    #             json.dump(data, output)
    
# wildcard_constraints:
#     sample='^.*(?=/)'


################################################################################
# Split by SSS barcodes
################################################################################


checkpoint extract_split_sss:
    wildcard_constraints:
        sample='|'.join(samples)
    input:
        # r1 = out_dir + "fastqs_split/{sample}.R1.fastq.gz",
        # r2 = out_dir + "fastqs_split/{sample}.R2.fastq.gz"
        r1 = out_dir + "fastqs_split/{sample}/{sample}.R1.{sfq}.fastq.gz",
        r2 = out_dir + "fastqs_split/{sample}/{sample}.R2.{sfq}.fastq.gz"
    output:
        temp(directory(temp_dir + 'split_SSS_tmp/{sample}/{sfq}')) # does not work if this is a temp dir
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stdout = out_dir + "logs/{sample}_{sfq}_sss.out",
        stderr = out_dir + "logs/{sample}_{sfq}_sss.err"
    params:
        # output_dir=out_dir + 'split_SSS/{sample}/',
        barcodes=lambda wildcards: config['samples'][wildcards.sample],
        sss_m=config['sss_mismatches'],
        rt_m=config['rt_mismatches'],
        rt_p=config['rt_primer']
    shell:
        '''
        python workflow/scripts/order_and_split_read_pairs.py \
        -r1 {input.r1} \
        -r2 {input.r2} \
        -o {output} \
        -b {params.barcodes} \
        -sss {params.sss_m} \
        -rt {params.rt_m} \
        -p {params.rt_p} > {log.stdout} 2> {log.stderr}

        for i in {output}/*.fastq; do
            pigz $i
        done 
        '''


# def extract_split_sss_output(wildcards):

#     outputs_i = glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards).output}/*/")
#     # print("ot", outputs_i)
#     outputs_i = [output.split('/')[-2] for output in outputs_i]
#     # print("outputs_i", outputs_i)
#     split_files = []
#     for i in outputs_i:
#         split_files.extend(expand([out_dir + 'split_SSS/{{sample}}/{{i}}.R1.ordered.fastq',
#                                    out_dir + 'split_SSS/{{sample}}/{{i}}.R1.ordered.fastq'],
#                                    sample=wildcards.sample))
#     # print("ess_split_file", split_files)
#     return split_files



# rule extract_split_sss_check:
#     input:
#         extract_split_sss_output
#     output:
#         touch(out_dir + 'logs/{sample}_sss_all_done.txt')



################################################################################
# Merge fastqs after SSS split finished
################################################################################


def merge_sss_fastqs_input(wildcards, r=1):

    outputs_i = glob.glob(f"{checkpoints.split_pe.get(**wildcards).output}/*R1*.fastq.gz")
    outputs_i = [output.split('/')[-1].split('.')[2] for output in outputs_i]

    split_files = []

    if r==1:

        split_files.extend(expand(temp_dir + f"split_SSS_tmp/{{sample}}/{{i}}/{{sss}}.R1.ordered.fastq.gz",
                                sample=wildcards.sample,
                                sss=wildcards.sss,
                                i=outputs_i))

    if r==2:

        split_files.extend(expand(temp_dir + f"split_SSS_tmp/{{sample}}/{{i}}/{{sss}}.R2.ordered.fastq.gz",
                                sample=wildcards.sample,
                                sss=wildcards.sss,
                                i=outputs_i))

    # print("summary:", split_files)
    return split_files



rule merge_sss_fastqs:
    input:
        r1 = lambda wc: merge_sss_fastqs_input(wc, r=1),
        r2 = lambda wc: merge_sss_fastqs_input(wc, r=2)
    output:
        r1 = temp_dir + "split_SSS/{sample}/{sss}.R1.ordered.fastq.gz",
        r2 = temp_dir + "split_SSS/{sample}/{sss}.R2.ordered.fastq.gz"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stderr = out_dir + "logs/{sample}/{sss}_merge_sss_fastqs.err"
    shell:
        '''
        cat {input.r1} > {output.r1} 2> {log.stderr}
        cat {input.r2} > {output.r2} 2>> {log.stderr}
        '''


