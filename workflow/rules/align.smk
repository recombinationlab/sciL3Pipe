
def align_input(wildcards, TRM_FILES=TRM_FILES):
    
    if trimmed_start:
        out = TRM_FILES.get(wildcards.sample).get(wildcards.sss)
    else:
        out = [out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
               out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.attached.fastq.gz"]
    return out



rule bwa_align_pe:
    '''
    -M 	Mark shorter split hits as secondary
    '''
    input:
        # fq_1=out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
        # fq_2=out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.attached.fastq.gz"
        align_input
    output:
        temp(out_dir + f"bwa_alignment/{{sample}}/{{sss}}.PE.bwa.{assembly}.bam")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.PE.bwa.{assembly}.log"
    params:
        rg="\'@RG\\tID:{sss}\\tSM:{sss}\\tPL:ILLUMINA\'"
    threads: 
        10
    shell:
        '''
        (bwa mem \
        -t {threads} \
        -M \
        -R {params.rg} \
        {bwa_index} \
        {input} \
        | \
        samtools view -b -F 256 - > {output}) &> {log}
        '''


rule bowtie2_align_pe:
    '''
    MapQ filter 20, -F 4 only mapped reads, -F 256 remove not primary alignment reads
    -F: Do not output alignments with any bits set in INT present in the FLAG field
    --maxins Bowtie2 by default has max pair insert size of 500
    '''
    input:
        # fq_1=out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
        # fq_2=out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.attached.fastq.gz"
        align_input
    output:
        temp(out_dir + f"bowtie2_alignment/{{sample}}/{{sss}}.PE.bowtie2.{assembly}.bam")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.PE.bowtie2.{assembly}.log"
    threads: 
        10
    shell:
        '''
        (bowtie2 \
        -p {threads} \
        -t \
        --maxins 2000 \
        --phred33 \
        -x {bowtie2_index} \
        -1 {input.fq_1} -2 {input.fq_2} | \
        samtools view -@ {threads} -b -F 256 - > {output}) &> {log}
        '''

rule sort_index_markdup_pe:
    '''
    should be collated coming from aligner
    fill in mate coordinates and insert size fields
    sort
    mark duplicates (performed in downstream pipeline)
    index
    '''
    input:
        out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.bam"
    output:
        fixmate=temp(out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.fixmate.bam"),
        bam=out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam",
        # mrkdup=out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam",
        idx=out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam.bai",
        flagstat=out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam.flagstat",
        idxstats=out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam.idxstats",
        stats=out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam.stats"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.sort_index_markdup_pe.err"
    threads:
        5
    shell:
        '''
        samtools fixmate -@ {threads} -m {input} {output.fixmate} 2> {log}

        samtools sort -@ {threads} -o {output.bam} {output.fixmate} 2>> {log}
        
        samtools index {output.bam} 2>> {log}

        samtools flagstat {output.bam} > {output.flagstat} 2>> {log}
        samtools idxstats {output.bam} > {output.idxstats} 2>> {log}
        samtools stats {output.bam} > {output.stats} 2>> {log}
        '''

rule qc_plots:
    '''
    Need to sort by name before running python script
    Will produce an html with qc plots, does not do any actual filtering, but plots should be used to filtering in downstream steps.
    '''
    input:
        out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam"
    output:
        report(out_dir + f"{aligner}_alignment/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.srt.bam.html", caption="../report/Alignment_QC.rst")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.{aligner}.{assembly}.qc_plot.log"
    params:
        hr='--hybrid_reference' if config['hybrid_reference'] else ''
    threads:
        10
    shell:
        '''
        python workflow/scripts/filter_bam.py \
        --edit_max 0 \
        --edit_min 0 \
        -i {input} \
        --insert_max 0 \
        --insert_min 0 \
        --mapq_max 255 \
        --mapq_min 20 \
        {params.hr} \
        --threads {threads} \
        --paired &> {log}
        '''


def align_output(wildcards):

    if trimmed_start:
        out_qc = list()
        for s in TRM_FILES.keys():
            if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                sss_keys = list(TRM_FILES.get(s).keys())
                for i in sss_keys:
                    out_qc.extend(expand(out_dir + f"{aligner}_alignment/{s}/{i}.PE.{aligner}.{assembly}.srt.bam.html"))
    else:
        output_i = glob.glob(f"{checkpoints.split_pe.get(**wildcards).output}/*R1*.fastq.gz")
        outputs_i = [output.split('/')[-1].split('.')[2] for output in output_i]

        # print("outputs_i", outputs_i)
        out_qc = []
        for i_wc in outputs_i:

            # print("output", glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards, sfq=i_wc).output}/*R[12].ordered.*"))
            output_j = glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards, sfq=i_wc).output}/*R[12].ordered.*")
            outputs_j = [output.split('/')[-1].split('.')[0] for output in output_j]
            # print(output_j)
            # print("outputs_j", outputs_j)
            out_qc.extend(expand(out_dir + f"{aligner}_alignment/{{sample}}/{{j}}.PE.{aligner}.{assembly}.srt.bam.html",
                                sample=wildcards.sample,
                                j=outputs_j))
        # print('out_qc', out_qc)
    return out_qc




# def qc_plots_output(wildcards):
#     checkpoint_output = checkpoints.extract_split_sss.get(**wildcards).output[0]
#     out = expand(out_dir + f"{aligner}_alignment/{{sample}}/{{i}}.PE.{aligner}.{assembly}.srt.bam.html",
#                sample=wildcards.sample,
#                i=glob_wildcards(os.path.join(checkpoint_output, "{i}.R1.ordered.fastq")).i)
#     if len(out) == 0:
#         out = expand(out_dir + f"{aligner}_alignment/{{sample}}/{{i}}.PE.{aligner}.{assembly}.srt.bam.html",
#                sample=wildcards.sample,
#                i=glob_wildcards(os.path.join(checkpoint_output, "{i}.R1.ordered.fastq.gz")).i)
#     # print('QC plots output:', out)
#     return out


checkpoint collect_and_generate_SSS_json:
    '''
    create SSS json used as input into split_bam

    merge_sss_bam_input does not work if this is not a checkpoint
    '''
    input:
        align_output
    output:
        out_dir + f"sss_json_{aligner}_{assembly}/{{sample}}/sss_samples.json"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stdout = out_dir + f"logs/{{sample}}.collect_and_generate_SSS_json_{aligner}_{assembly}.log",
        stderr = out_dir + f"logs/{{sample}}.collect_and_generate_SSS_json_{aligner}_{assembly}.err"
    params:
        bam_dir = out_dir + f"{aligner}_alignment/{{sample}}",
        ends = f"{aligner}.{assembly}.srt.bam"
    shell:
        '''
        python sss_bam2json.py --bam_dir {params.bam_dir} --out_file {output} --ends {params.ends} > {log.stdout} 2> {log.stderr}
        '''
