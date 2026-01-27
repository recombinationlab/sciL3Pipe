
def merge_sss_bam_input(wildcards, SSS_FILES=SSS_FILES):
    
    if sss_bam_start:
        out = SSS_FILES.get(wildcards.sample).get(wildcards.sss)
    elif trimmed_start:
        out = out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam"
    else:
        paths = []
        # print(rules.collect_and_generate_SSS_json.output)
        for i in expand(checkpoints.collect_and_generate_SSS_json.get(**wildcards).output[0]):
            # expand(out_dir + f"sss_json/{wildcards.sample}/sss_samples.json"):
            SSS_FILES = json.load(open(i))
            # print(SSS_FILES.get(wildcards.sample).get(wildcards.sss))
            paths.extend(SSS_FILES.get(wildcards.sample).get(wildcards.sss))
        out = []
        # ensure all paths end with .srt.bam
        for bam in paths:
            if bam.endswith('.srt.bam'):
                out.append(bam)

    # print('out test', out)
    return out


rule merge_sss_bam:
    '''
    Combine BAM files with the same SSS, for example if sequenced multiple times.
    Avoids the need to run previous steps and files can be combined at this stage from multiple folders.
    '''
    input:
        merge_sss_bam_input
        # lambda wildcards: SSS_FILES.get(wildcards.sample).get(wildcards.sss)
    output:
        out_dir + f"merged_sss_bam/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + "logs/{sample}/{sss}.merge_sss_bam.err"
    params:
        dir_out = out_dir + f"merged_sss_bam/{{sample}}/"
    threads:
        10
    shell:
        '''
        count=$(echo '{input}' | awk -F' ' '{{print NF}}')
        mkdir -p {params.dir_out} 2> {log}

        if [[ $count -gt 1 ]]
        then
            samtools merge --threads {threads} - {input} | samtools sort --threads {threads} -o {output} 2>> {log}
            samtools index {output} 2>> {log}
        else
            ln -sr {input} {output} 2>> {log}
        fi
        '''
 

checkpoint split_bam:
    '''
    Split bam file into individual single cell bam files based on barcode
    '''
    wildcard_constraints:
        sample='|'.join(samples),
        sss='|'.join(sss_all) if len(sss_all) > 0 else r'[^/]+'
    input:
        out_dir + f"merged_sss_bam/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam" #\
        # if sss_bam_start else \
        # out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam"        
    output:
        directory(out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.{aligner}.{assembly}.split_bam.log"
    params: 
        tn5_b=config['tn5_barcodes'],
        l_b=config['ligation_barcodes'],
        cutoff=config['cutoff']
    shell:
        '''
        python workflow/scripts/split_bam.py \
        -i {input} \
        -o {output} \
        -tb {params.tn5_b} \
        -lb {params.l_b} \
        --cutoff {params.cutoff} &> {log}
        '''

# TODO: add split_bam png to report 

rule index_bam:
    input:
        out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam"
    output:
        out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam.bai"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stdout = out_dir + f"logs/{{sample}}/{{sss}}/{{i}}.{aligner}.{assembly}.index_bam.log",
        stderr = out_dir + f"logs/{{sample}}/{{sss}}/{{i}}.{aligner}.{assembly}.index_bam.err"
    threads:
        1
    group:
        "groupjob_idx"
    shell:
        '''
        samtools index -@ {threads} {input} {output} > {log.stdout} 2> {log.stderr}
        '''

