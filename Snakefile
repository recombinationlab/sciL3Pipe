
import datetime
from pathlib import Path
import yaml
from shutil import which

# SAMPLES = ["yi293_GAACCG"]
# TODO: remove markdup and implement umi_tools
# TODO: update filter_bam.py rule (don't need to collate sort anymore)

v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')

INFO = ""

try:
    out_dir = config['output_dir']
    INFO += 'All data will be written to: ' + out_dir + '\n'
except:
    out_dir = ''
    INFO += 'Defaulting to working directory as output directory' + '\n'

try:
    aligner = config['aligner']
except:
    sys.exit('ERROR: Need to specify aligner')

if aligner == 'bwa':
    try:
        bwa_index = config['bwa_index'][config['assembly']]
    except:
        sys.exit('ERROR: BWA index not specified in config.yaml')
elif aligner == 'bowtie2':
    try:
        bowtie2_index = config['bowtie2_index'][config['assembly']]
    except:
        sys.exit('ERROR: Bowtie2 index not specified in config.yaml')

try:
    assembly = config['assembly']
except:
    sys.exit('ERROR: Need to specify assembly e.g. "hg38" or "mm10"')


################################################################################
# Load samples 
################################################################################

# Fastq files
try:
    FILES = json.load(open(config['fastqs']))
except:
    FILES = None

# Samples should start with a general name e.g. yi292, followed by
# SSS file containing barcodes
ALL_SAMPLES = []
POST_SSS_SAMPLES = []
for samp in config['samples']:
    ALL_SAMPLES.append(samp)
    with open(config['samples'][samp]) as sss:
        for line in sss:
            POST_SSS_SAMPLES.append(samp + '_' + line.rstrip())

INFO += f'Samples to process {ALL_SAMPLES} \n'
INFO += f'SSS split samples to process {POST_SSS_SAMPLES} \n'

################################################################################
# Out files
################################################################################

MERGE = expand([out_dir + "merged/{sample}.R1.fastq.gz",
                out_dir + "merged/{sample}.R2.fastq.gz"],
                 sample=ALL_SAMPLES)


TRIM = expand([out_dir + "trimmed/{sample}/{sss}.R1.trimmed.fq.gz",
               out_dir + "trimmed/{sample}/{sss}.R2.trimmed.fq.gz"],
                 sample=ALL_SAMPLES, sss=POST_SSS_SAMPLES)

ATTACHED = expand([out_dir + "trimmed/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
                   out_dir + "trimmed/{sample}/{sss}.R2.trimmed.attached.fastq.gz",
                   out_dir + "trimmed/{sample}/{sss}.noME.R1-2.fastq.gz"],
                   sample=ALL_SAMPLES, sss=POST_SSS_SAMPLES)

G_COV = expand(out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}.genome_coverage.txt",
               sample = ALL_SAMPLES, sss = POST_SSS_SAMPLES)

COLLISIONS = expand([out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}.collision.txt",
                     out_dir + f"genome_coverage_{aligner}/single_cells/{{sample}}/{{sss}}.collision.txt"],
                     sample= ALL_SAMPLES, sss = POST_SSS_SAMPLES)

SUMMARY = [out_dir + f"genome_coverage_{aligner}/single_cells/single_cell_summary_all.txt"]

QC = expand([out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.collate.bam.NM.png",
             out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.collate.bam.MAPQ.png"],
             sample = ALL_SAMPLES, sss = POST_SSS_SAMPLES)

################################################################################
################################################################################
# Execute before workflow starts
################################################################################
################################################################################
onstart:
    print(INFO)
################################################################################
################################################################################
# Rule all
################################################################################
################################################################################
if FILES == None:
    rule all:
        input: 
            G_COV + SUMMARY + COLLISIONS + QC
else:
    rule all:
        input: 
            MERGE + TRIM + ATTACHED + G_COV + SUMMARY + COLLISIONS + QC


################################################################################
# Log
################################################################################

rule log_config:
    '''Copy config.yaml and place in logs folder with the date run
    '''
    output:
        out_dir + "logs/config_" + run_date + "yaml"
    run:
        with open(output[0], 'w') as out:
            yaml.dump(config, out, default_flow_style=False)

################################################################################
# Rename merge
################################################################################

rule merge_fastqs_pe:
    # wildcard_constraints:
    #     sample='^.*(?=/)'
    input: 
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: 
        r1 = out_dir + "merged/{sample}.R1.fastq.gz",
        r2 = out_dir + "merged/{sample}.R2.fastq.gz"
    conda:
        "envs/sciStrand_env.yaml"
    shell:
        ''' 
        count_1=$(echo '{input.r1}' | awk -F' ' '{{print NF}}')
        
        if [[ $count_1 -gt 1 ]]
        then
            cat {input.r1} > {output.r1}
        else
            ln -sr {input.r1} {output.r1}
        fi

        count_2=$(echo '{input.r2}' | awk -F' ' '{{print NF}}')
        
        if [[ $count_2 -gt 1 ]]
        then
            cat {input.r2} > {output.r2}
        else
            ln -sr {input.r2} {output.r2}
        fi
        '''

################################################################################
# Split by SSS barcodes
################################################################################

rule extract_split_sss:
    input:
        r1 = out_dir + "merged/{sample}.R1.fastq.gz",
        r2 = out_dir + "merged/{sample}.R2.fastq.gz"
    output:
        expand([out_dir + 'split_SSS/{{sample}}/{sss}.R1.ordered.fastq',
                out_dir + 'split_SSS/{{sample}}/{sss}.R2.ordered.fastq'],
                sss=POST_SSS_SAMPLES)
    conda:
        "envs/sciStrand_env.yaml"
    log:
        out_dir + "logs/{sample}_sss.log"
    params:
        output_dir=out_dir + 'split_SSS/{sample}/',
        barcodes=lambda wildcards: config['samples'][wildcards.sample],
        sss_m=config['sss_mismatches'],
        rt_m=config['rt_mismatches'],
        rt_p=config['rt_primer']
    shell:
        '''
        python scripts/order_and_split_read_pairs.py \
        -r1 {input.r1} \
        -r2 {input.r2} \
        -o {params.output_dir} \
        -b {params.barcodes} \
        -sss {params.sss_m} \
        -rt {params.rt_m} \
        -p {params.rt_p} &> {log}
        '''


rule bgzip_sss_fastq:
    '''
    bgzip output of extract_split_sss checkpoint
    pigz - Allow up to n compression threads (default is the
    number of online processors, or 8 if unknown)
    '''
    input:
        out_dir + 'split_SSS/{sample}/{sss}.R1.ordered.fastq',
        out_dir + 'split_SSS/{sample}/{sss}.R2.ordered.fastq'
    output:
        out_dir + 'split_SSS/{sample}/{sss}.R1.ordered.fastq.gz',
        out_dir + 'split_SSS/{sample}/{sss}.R2.ordered.fastq.gz'
    conda:
        "envs/sciStrand_env.yaml"
    threads:
        10
    shell: 
        '''
        pigz -p {threads} {output}
        '''

################################################################################
# Trim barcodes
################################################################################

rule cutadapt:
    '''
    Trim Tn5 sequences
    The --info-file, --rest-file and --wildcard-file options write out information only from the first read.
    '''
    input:
        [out_dir + "split_SSS/{sample}/{sss}.R1.ordered.fastq.gz",
         out_dir + "split_SSS/{sample}/{sss}.R2.ordered.fastq.gz"]
    output:
        fastq1=out_dir + "trimmed/{sample}/{sss}.R1.trimmed.fq.gz",
        fastq2=out_dir + "trimmed/{sample}/{sss}.R2.trimmed.fq.gz",
        info=out_dir + "trimmed/{sample}/{sss}.R1.info.tab",
        qc=out_dir + "trimmed/{sample}/{sss}.trim.qc.txt"
    log:
        out_dir + "logs/{sample}/{sss}_cutadapt_r1.log"
    conda:
        "envs/sciStrand_env.yaml"
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
        out_dir + "trimmed/{sample}/{sss}.R1.info.tab"
    output:
        out_dir + "trimmed/{sample}/{sss}.R1.bc1.bc2.txt.gz"
    conda:
        "envs/sciStrand_env.yaml"
    shell:
        '''
        awk -F $'\t' '{{if($2>=0) {{print substr($5,length($5)-20,7), substr($5,length($5)-7,8);}} else print "NA","NA"}}' \
        {input} | gzip - > {output}
        '''


rule attach_tn5_barcodes:
    input:
        fastq1=out_dir + "trimmed/{sample}/{sss}.R1.trimmed.fq.gz",
        fastq2=out_dir + "trimmed/{sample}/{sss}.R2.trimmed.fq.gz",
        bc1_bc2=out_dir + "trimmed/{sample}/{sss}.R1.bc1.bc2.txt.gz"
    output:
        [out_dir + "trimmed/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
        out_dir + "trimmed/{sample}/{sss}.R2.trimmed.attached.fastq.gz",
        out_dir + "trimmed/{sample}/{sss}.noME.R1-2.fastq.gz"]
    conda:
        "envs/sciStrand_env.yaml"
    log:
        out_dir + "logs/{sample}/{sss}.tn5_ligation_barcodes.log"
    params:
        dir=out_dir + "trimmed/{sample}",
        to_gzip=[out_dir + "trimmed/{sample}/{sss}.R1.trimmed.attached.fastq",
                 out_dir + "trimmed/{sample}/{sss}.R2.trimmed.attached.fastq",
                 out_dir + "trimmed/{sample}/{sss}.noME.R1-2.fastq"],
        tn5_b=config['tn5_barcodes'],
        l_b=config['ligation_barcodes'],
        tn5_m=config['tn5_mismatches'],
        l_m=config['ligation_mismatches']
    threads:
        10
    shell:
        '''
        python scripts/parse_tn5_ligation_barcodes.py \
        -r1 {input.fastq1} \
        -r2 {input.fastq2} \
        -o {params.dir} \
        -bc {input.bc1_bc2} \
        -bt {params.tn5_b} \
        -bl {params.l_b} \
        -tn5 {params.tn5_m} \
        -l {params.l_m} &> {log}

        pigz {params.to_gzip}
        '''

################################################################################
# Align rule
################################################################################

rule bwa_align_pe:
    '''
    -M 	Mark shorter split hits as secondary
    '''

    input:
        fq_1=out_dir + "trimmed/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
        fq_2=out_dir + "trimmed/{sample}/{sss}.R2.trimmed.attached.fastq.gz"
    output:
        out_dir + f"alignments/{{sample}}/{{sss}}.PE.bwa.{assembly}.bam"
    threads: 10
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.PE.bwa.{assembly}.log"
    conda:
        "envs/sciStrand_env.yaml"
    params:
        rg="\'@RG\\tID:{sss}\\tSM:{sss}\\tPL:ILLUMINA\'"
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
        fq_1=out_dir + "trimmed/{sample}/{sss}.R1.trimmed.attached.fastq.gz",
        fq_2=out_dir + "trimmed/{sample}/{sss}.R2.trimmed.attached.fastq.gz"
    output:
        out_dir + f"alignments/{{sample}}/{{sss}}.PE.bowtie2.{assembly}.bam"
    threads: 10
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.PE.bowtie2.{assembly}.log"
    conda:
        "envs/sciStrand_env.yaml"
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
    mark duplicates
    index
    '''
    input:
        out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.bam"
    output:
        fixmate=temp(out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.fixmate.bam"),
        bam=temp(out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.sorted.bam"),
        mrkdup=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam",
        idx=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam.bai",
        flagstat=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam.flagstat",
        idxstats=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam.idxstats",
        stats=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam.stats"
    conda:
        "envs/sciStrand_env.yaml"
    threads:
        10
    shell:
        '''
        samtools fixmate -@ {threads} -m {input} {output.fixmate}

        samtools sort -@ {threads} -o {output.bam} {output.fixmate}
        
        samtools markdup -@ {threads} {output.bam} {output.mrkdup}
        samtools index {output.mrkdup}

        samtools flagstat {output.mrkdup} > {output.flagstat}
        samtools idxstats {output.mrkdup} > {output.idxstats}
        samtools stats {output.mrkdup} > {output.stats}
        '''

rule qc_plots:
    '''
    Need to sort by name before running python script
    Will produce 4 png plots
    '''
    input:
        out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam"
    output:
        collate=temp(out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.collate.bam"),
        png1=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.collate.bam.NM.png",
        png2=out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.collate.bam.MAPQ.png"
    conda:
        "envs/sciStrand_env.yaml"
    params:
        hr='--hybrid_reference' if config['hybrid_reference'] else ''
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.{aligner}.qc_plot.log"
    threads:
        10
    shell:
        '''
        samtools collate -@ {threads} -o {output.collate} {input}

        python scripts/filter_bam.py \
        --edit_max 0 \
        --edit_min 0 \
        -i {output.collate} \
        --insert_max 0 \
        --insert_min 0 \
        --mapq_max 255 \
        --mapq_min 20 \
        {params.hr} \
        --paired &> {log}
        '''

################################################################################
# Split BAM into single cell BAMs
################################################################################

checkpoint split_bam:
    '''
    Split bam file into individual single cell bam files based on barcode
    '''
    input:
        out_dir + f"alignments/{{sample}}/{{sss}}.PE.{aligner}.{assembly}.markdup.bam"
    output:
        directory(out_dir + f"split_bam_{aligner}/{{sample}}/{{sss}}")
    conda:
        "envs/sciStrand_env.yaml"
    log:
        out_dir + f"logs/{{sample}}/{{sss}}.{aligner}.split_bam.log"
    params: 
        tn5_b=config['tn5_barcodes'],
        l_b=config['ligation_barcodes'],
        cutoff=config['cutoff']
    shell:
        '''
        python scripts/split_bam.py \
        -i {input} \
        -o {output} \
        -tb {params.tn5_b} \
        -lb {params.l_b} \
        --cutoff {params.cutoff} &> {log}
        '''


################################################################################
# Summary metrics
################################################################################

rule bed_files:
    input:
        out_dir + f"split_bam_{aligner}/{{sample}}/{{sss}}/{{i}}.bam"
    output:
        idx=out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}/{{i}}.bam.idx",
        bed=out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.bed",
        bed_R1=out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.bed",
        bed_R1_human=out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.human.bed",
        bed_q=out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.human.Qgt0.bed",
        bed_q_unq=out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.human.Qgt0.uniq.bed",
        uniq1=temp(out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.human.Qgt0.uniq1.bed"),
        uniq2=temp(out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.human.Qgt0.uniq2.bed")
    conda:
        "envs/sciStrand_env.yaml"
    shell:
        '''
        samtools index {input}
        samtools idxstats {input} > {output.idx}

        bedtools bamtobed -i {input} > {output.bed}
        # Get read 1
        grep "/1" {output.bed} > {output.bed_R1}
        # Get human (without chr) counts
        awk '{{if ($1 !~ /^chr/) print}}' {output.bed_R1} > {output.bed_R1_human}
        awk '{{if ($5 > 0) print}}' {output.bed_R1_human} > {output.bed_q}
        
        awk '!seen[$1, $2]++' {output.bed_q} > {output.uniq1}
        awk '!seen[$1, $3]++' {output.bed_q} > {output.uniq2}
        bedtools intersect -a {output.uniq1} -b {output.uniq2} -f 1 -F 1 > {output.bed_q_unq}
        '''
   
def split_bam_output(wildcards):
    checkpoint_output = checkpoints.split_bam.get(**wildcards).output[0]
    return expand(out_dir + f"split_bam_{aligner}/{{sample}}/{{sss}}/{{i}}.bam",
               sample=wildcards.sample,
               sss=wildcards.sss,
               i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bam")).i)

def bed_gcov(wildcards):
    checkpoint_output = checkpoints.split_bam.get(**wildcards).output[0]
    return expand(out_dir + f"bed_files_{aligner}/{{sample}}/{{sss}}/{{i}}.R1.human.Qgt0.bed",
               sample=wildcards.sample,
               sss=wildcards.sss,
               i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bam")).i)

def split_bam_idx(wildcards):
    checkpoint_output = checkpoints.split_bam.get(**wildcards).output[0]
    return expand(out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}/{{i}}.bam.idx",
           sample=wildcards.sample,
           sss=wildcards.sss,
           i=glob_wildcards(os.path.join(checkpoint_output, "{i}.bam")).i)


rule genome_coverage:
    '''
    Requires bai and idx created in rule bed_files
    '''
    input:
        bam=split_bam_output,
        bed_q=bed_gcov,
        idx=split_bam_idx
    output:
        out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}.genome_coverage.txt"
    conda:
        "envs/sciStrand_env.yaml"
    shell:
        '''
        for i in {input.bam}; do

            bn=`basename ${{i%.*}}`
            echo $bn >> {output}
            # Get human (without chr) genome coverage
            awk '{{if ($1 !~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}/{wildcards.sample}/{wildcards.sss}/${{bn}}.bam.idx \
            | awk '{{sum += $3}} END {{print sum}}' >> {output}
                                    
            wc -l < {out_dir}bed_files_{aligner}/{wildcards.sample}/{wildcards.sss}/${{bn}}.R1.human.Qgt0.bed >> {output}
            wc -l < {out_dir}bed_files_{aligner}/{wildcards.sample}/{wildcards.sss}/${{bn}}.R1.human.Qgt0.uniq.bed >> {output}
        done
        '''



rule summary:
    input:
        expand(out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}.genome_coverage.txt", 
               sample=ALL_SAMPLES, sss=POST_SSS_SAMPLES)
    output:
        out_dir + "genome_coverage_{aligner}/single_cells/single_cell_summary_all.txt"
    conda:
        "envs/sciStrand_env.yaml"
    shell:
        '''
        for i in {input}; do
            awk 'NR%4{{printf "%s ",$0;next;}}1' $i | awk '{{print $1, $2, $3, $4, $3/$4}}' >> {output}
        done
        '''

    
rule collisions:
    '''
    Requires bai and idx created in rule bed_files
    '''
    input:
        bam=split_bam_output,
        idx=split_bam_idx
    output:
        o1=out_dir + f"genome_coverage_{aligner}/{{sample}}/{{sss}}.collision.txt",
        o2=out_dir + f"genome_coverage_{aligner}/single_cells/{{sample}}/{{sss}}.collision.txt"
    shell:
        '''
        for i in {input.bam}; do

            bn=`basename $i`
            echo ${{bn%.*}} >> {output.o1}
            awk '{{if ($1 !~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}/{wildcards.sample}/{wildcards.sss}/${{bn}}.idx \
            | awk '{{sum += $3}} END {{print sum}}' >> {output.o1}
            awk '{{if ($1 ~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}/{wildcards.sample}/{wildcards.sss}/${{bn}}.idx \
            | awk '{{sum += $3}} END {{print sum}}' >> {output.o1}
            
            awk 'NR%3{{printf "%s ",$0;next;}}1' {output.o1} | awk '{{print $1, $2, $3, ($2+$3), $2/($2+$3), $3/($2+$3)}}' > {output.o2}

            cat {output.o2} >> {out_dir}genome_coverage_{aligner}/single_cells/collision.txt
        done
        '''

################################################################################
################################################################################
# After succesful run
################################################################################
################################################################################
# Delete all temporary directories in case temp() did not remove them
onsuccess:
    assert which('multiqc') is not None, 'ERROR: multiqc not found'
    #run multiqc after pipeline completes
    multiqc_out = Path(out_dir + "multiQC/multiqc_report.html")
    if not multiqc_out.exists():
        shell('multiqc {out_dir} -o {out_dir}multiQC')