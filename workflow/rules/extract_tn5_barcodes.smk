
rule cutadapt:
    '''
    Trim Tn5 sequences
    The --info-file, --rest-file and --wildcard-file options write out information only from the first read.
    '''
    input:
        [out_dir + "split_SSS/{sample}/{sss}.R1.ordered.fastq.gz",
         out_dir + "split_SSS/{sample}/{sss}.R2.ordered.fastq.gz"]
    output:
        fastq1=temp(out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.fq.gz"),
        fastq2=temp(out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.fq.gz"),
        info=temp(out_dir + "trimmed_fastq/{sample}/{sss}.R1.info.tab"),
        qc=out_dir + "trimmed_fastq/{sample}/{sss}.trim.qc.txt"
    log:
        out_dir + "logs/{sample}/{sss}_cutadapt_r1.log"
    conda:
        "../envs/sciL3_env.yaml"
    threads: 
        10
    shell:
        '''
        cutadapt \
        -g "AGATGTGTATAAGAGACAG;e=0.2;o=19" \
        -G "AGATGTGTATAAGAGACAG;e=0.2;o=13" \
        -A "CTGTCTCTTATACACATCT;e=0.2;o=13" \
        --info-file {output.info} \
        -o {output.fastq1} \
        -p {output.fastq2} \
        -j {threads} \
        {input} 1> {output.qc} 2> {log}
        '''

rule extract_tn5_barcodes:
    '''
    Extract Tn5 barcodes from info file
    '''
    input:
        out_dir + "trimmed_fastq/{sample}/{sss}.R1.info.tab"
    output:
        temp(out_dir + "trimmed_fastq/{sample}/{sss}.R1.bc1.bc2.txt.gz")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stderr = out_dir + "logs/{sample}/{sss}.extract_tn5_barcodes.err"
    shell:
        '''
        awk -F $'\t' '{{if($2>=0) {{print substr($5,length($5)-20,7), substr($5,length($5)-7,8);}} else print "NA","NA"}}' \
        {input} | gzip - > {output} 2> {log.stderr}
        '''

#TODO add sequencing run date to the name of these trimmed fqs, also add it to the split_SSS fq files
rule attach_tn5_barcodes:
    input:
        fastq1=out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.fq.gz",
        fastq2=out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.fq.gz",
        bc1_bc2=out_dir + "trimmed_fastq/{sample}/{sss}.R1.bc1.bc2.txt.gz"
    output:
        [out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
        out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.attached.fastq.gz",
        temp(out_dir + "trimmed_fastq/{sample}/{sss}.noME.R1-2.fastq.gz"),
        out_dir + "logs/{sample}/{sss}.tn5_ligation_barcodes.log"]
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stdout = out_dir + "logs/{sample}/{sss}.tn5_ligation_barcodes.log",
        stderr = out_dir + "logs/{sample}/{sss}.tn5_ligation_barcodes.err"
    params:
        dir=out_dir + "trimmed_fastq/{sample}",
        to_gzip=[out_dir + "trimmed_fastq/{sample}/{sss}.R1.trimmed.attached.fastq",
                 out_dir + "trimmed_fastq/{sample}/{sss}.R2.trimmed.attached.fastq",
                 out_dir + "trimmed_fastq/{sample}/{sss}.noME.R1-2.fastq"],
        tn5_b=config['tn5_barcodes'],
        l_b=config['ligation_barcodes'],
        tn5_m=config['tn5_mismatches'],
        l_m=config['ligation_mismatches']
    shell:
        '''
        python workflow/scripts/parse_tn5_ligation_barcodes.py \
        -r1 {input.fastq1} \
        -r2 {input.fastq2} \
        -o {params.dir} \
        -bc {input.bc1_bc2} \
        -bt {params.tn5_b} \
        -bl {params.l_b} \
        -tn5 {params.tn5_m} \
        -l {params.l_m} > {log.stdout} 2> {log.stderr}

        pigz {params.to_gzip} 2>> {log.stderr}
        '''


# rule summarise_tn5_ligation_barcodes_log:
#     input:
#         out_dir + "logs/{sample}/{sss}.tn5_ligation_barcodes.log"
#     output:
