
rule bed_files:
    '''
    Get read counts (R1) and unique read counts by star and end position of R1
    '''
    input:
        bam=out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam",
        bai=out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam.bai"
    output:
        idx=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam.idx",
        bed=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bed"),
        bed_R1=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.bed"),
        bed_R1_human=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.bed"),
        bed_q=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.Qgt0.bed"),
        bed_q_unq=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.Qgt0.uniq.bed"),
        uniq1=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.Qgt0.uniq1.bed"),
        uniq2=temp(out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.Qgt0.uniq2.bed")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + "logs/{sample}/{sss}/{i}.bed_files.log"
    threads:
        1
    # container:
    #     None
    group:
        "groupjob_bed"
    params:
        sf = summary_for
    shell:
        '''
        samtools idxstats {input.bam} > {output.idx} 2> {log}

        bedtools bamtobed -i {input.bam} > {output.bed} 2>> {log}
        # Get read 1
        grep "/1" {output.bed} > {output.bed_R1} 2>> {log}
        if [ {params.sf} == "wprefix" ]; then
            # Get human (without chr) counts
            awk '{{if ($1 ~ /^chr/) print}}' {output.bed_R1} > {output.bed_R1_human} 2>> {log}
            awk '{{if ($5 > 0) print}}' {output.bed_R1_human} > {output.bed_q} 2>> {log}
        elif [ {params.sf} == "woprefix" ]; then
            # Get human (without chr) counts
            awk '{{if ($1 !~ /^chr/) print}}' {output.bed_R1} > {output.bed_R1_human} 2>> {log}
            awk '{{if ($5 > 0) print}}' {output.bed_R1_human} > {output.bed_q} 2>> {log}
        fi
        # Unique by start and end poistion of R1
        awk '!seen[$1, $2]++' {output.bed_q} > {output.uniq1} 2>> {log}
        awk '!seen[$1, $3]++' {output.bed_q} > {output.uniq2} 2>> {log}
        bedtools intersect -a {output.uniq1} -b {output.uniq2} -f 1 -F 1 > {output.bed_q_unq} 2>> {log}
        '''
   


rule genome_coverage:
    '''
    Requires bai and idx created in rule bed_files
    '''
    input:
        bam=out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam",
        bed_q=out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.Qgt0.bed",
        bed_q_unq=out_dir + f"bed_files_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.R1.{summary_for}.Qgt0.uniq.bed",
        idx=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam.idx"
    output:
        out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.{summary_for}.genome_coverage.txt"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + "logs/{sample}/{sss}/{i}.genome_coverage.log"
    threads:
        1
    # container:
    #     None
    group:
        "groupjob_gc"
    params:
        sf = summary_for
    shell:
        '''
        for i in {input.bam}; do

            bn=`basename ${{i%.*}}`
            echo $bn >> {output}
            if [ {params.sf} == "wprefix" ]; then
                # Get mouse (with chr) genome coverage
                awk '{{if ($1 ~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}_{assembly}/{wildcards.sample}/{wildcards.sss}/${{bn}}.bam.idx \
                | awk '{{sum += $3}} END {{print sum}}' >> {output} 2>> {log}
            elif [ {params.sf} == "woprefix" ]; then
                # Get human (without chr) genome coverage
                awk '{{if ($1 !~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}_{assembly}/{wildcards.sample}/{wildcards.sss}/${{bn}}.bam.idx \
                | awk '{{sum += $3}} END {{print sum}}' >> {output} 2>> {log}
            fi
            wc -l < {out_dir}bed_files_{aligner}_{assembly}/{wildcards.sample}/{wildcards.sss}/${{bn}}.R1.{summary_for}.Qgt0.bed >> {output} 2>> {log}
            wc -l < {out_dir}bed_files_{aligner}_{assembly}/{wildcards.sample}/{wildcards.sss}/${{bn}}.R1.{summary_for}.Qgt0.uniq.bed >> {output} 2>> {log}
        done
        '''



    
rule collisions:
    '''
    Requires bai and idx created in rule bed_files
    '''
    input:
        bam=out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam",
        idx=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.bam.idx",
        gc=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.{summary_for}.genome_coverage.txt"
    output:
        o1=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.collision.txt",
        o2=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{{sss}}/{{i}}.single_cells.collision.txt"
    log:
        out_dir + "logs/{sample}/{sss}/{i}.collisions.log"
    threads:
        1
    # container:
    #     None
    group:
        "groupjob_col"
    shell:
        '''
        mkdir -p {out_dir}genome_coverage_{aligner}_{assembly}/{wildcards.sample}/single_cells
        for i in {input.bam}; do

            bn=`basename $i`
            echo ${{bn%.*}} >> {output.o1}
            awk '{{if ($1 !~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}_{assembly}/{wildcards.sample}/{wildcards.sss}/${{bn}}.idx \
            | awk '{{sum += $3}} END {{print sum}}' >> {output.o1} 2>> {log}
            awk '{{if ($1 ~ /^chr/) print}}' {out_dir}genome_coverage_{aligner}_{assembly}/{wildcards.sample}/{wildcards.sss}/${{bn}}.idx \
            | awk '{{sum += $3}} END {{print sum}}' >> {output.o1} 2>> {log}
            
            awk 'NR%3{{printf "%s ",$0;next;}}1' {output.o1} | awk '{{print $1, $2, $3, ($2+$3), $2/($2+$3), $3/($2+$3)}}' > {output.o2} 2>> {log}

            cat {output.o2} >> {out_dir}genome_coverage_{aligner}_{assembly}/{wildcards.sample}/single_cells/collision.txt 2>> {log}
        done
        '''




def summary_input(wildcards):
    if sss_bam_start:
        split_files = []
        for s in SSS_FILES.keys():
            if s in wildcards.sample: # this will prevent running all samples, rather run just the wildcard sample
                sss_keys = list(SSS_FILES.get(s).keys())
                for i in sss_keys:
                    # print('i', i)
                    # split_files.extend([out_dir + f"merged_sss_bam/{s}/{sss}.PE.{aligner}.{assembly}.srt.bam"])
                    output_j = glob.glob(f"{checkpoints.split_bam.get(sample=s, sss=i).output}/*.bam")
                    # print('output_j', output_j)
                    outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]
                    # print('outputs_j', outputs_j)
                    for j in outputs_j:
                        if config['generate_summary']:
                            split_files.extend(expand(out_dir + f"genome_coverage_{aligner}_{assembly}/{s}/{i}/{j}.{summary_for}.genome_coverage.txt"))
                        else:
                            split_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}/{s}/{i}/{j}.bam.bai"))
                        if config['run_collisions']:
                            split_files.extend(expand(out_dir + f"genome_coverage_{aligner}_{assembly}/{s}/{i}/{j}.single_cells.collision.txt"))
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
                        if config['generate_summary']:
                            split_files.extend(expand(out_dir + f"genome_coverage_{aligner}_{assembly}/{s}/{i}/{j}.{summary_for}.genome_coverage.txt"))
                        else:
                            split_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}/{s}/{i}/{j}.bam.bai"))
                        if config['run_collisions']:
                            split_files.extend(expand(out_dir + f"genome_coverage_{aligner}_{assembly}/{s}/{i}/{j}.single_cells.collision.txt"))

    else:

        outputs_i_wc = glob.glob(f"{checkpoints.split_pe.get(**wildcards).output}/*R1*.fastq.gz")
        outputs_i_wc = [output.split('/')[-1].split('.')[2] for output in outputs_i_wc]
        # print('outputs_i_wc', outputs_i_wc)

        split_files = []
        for i_wc in outputs_i_wc:
            # print('i_wc', i_wc)
            outputs_i = glob.glob(f"{checkpoints.extract_split_sss.get(**wildcards, sfq=i_wc).output}/*R1.ordered.*")
            outputs_i = [output.split('/')[-1].split('.')[0] for output in outputs_i]
            
            # print('outputs_i', outputs_i)
            for i in outputs_i:
                output_j = glob.glob(f"{checkpoints.split_bam.get(**wildcards, sss=i).output}/*.bam")
                outputs_j = [output.split('/')[-1].split('.bam')[0] for output in output_j]

                for j in outputs_j:
                    if config['generate_summary']:
                        split_files.extend(expand(out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{i}/{j}.{summary_for}.genome_coverage.txt",
                                                  sample=wildcards.sample))
                    else:
                        split_files.extend(expand(out_dir + f"split_bam_{aligner}_{assembly}/{{sample}}/{i}/{j}.bam.bai",
                                                  sample=wildcards.sample))
                    if config['run_collisions']:
                        split_files.extend(expand(out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/{i}/{j}.single_cells.collision.txt",
                                                  sample=wildcards.sample))
        # print("summary:", split_files)

    return split_files





rule summary_list:
    '''
    To avoid bash character limit, save all file names to file and use as input
    '''
    input:
        summary_input
    output:
        out_dir + f"files_{aligner}_{assembly}/{{sample}}_{summary_for}.list"
    run:
        with open(output[0], 'w') as out:
            out.write('\n'.join(input))


rule summary:
    input:
        out_dir + f"files_{aligner}_{assembly}/{{sample}}_{summary_for}.list"
    output:
        out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/single_cells/single_cell_summary_all_{summary_for}.txt"
    conda:
        "../envs/sciL3_env.yaml"
    log:
        out_dir + f"logs/{{sample}}.{aligner}_{assembly}.summary.log"
    params:
        create_dir=out_dir + f"genome_coverage_{aligner}_{assembly}/{{sample}}/single_cells"
    shell:
        '''
        mkdir -p {params.create_dir}
        if [ ! -s {input} ]; then
            touch {output}
        else
            for i in $(< {input}); do
                if [[ $i == *.genome_coverage.txt ]]; then
                    awk 'NR%4{{printf "%s ",$0;next;}}1' $i | awk '{{if($4==0) {{print $1, $2, $3, $4, 0}} else {{print $1, $2, $3, $4, $3/$4}}}}' >> {output} 2>> {log}
                fi
            done
        fi
        '''



def sss_log_input():
    result = []

    if SSS_FILES == None:
        samples = ALL_SAMPLES
    else:
        samples = SSS_SAMPLES
        
    for s in samples:
        result.append([out_dir + f"files_{aligner}_{assembly}/{s}_{summary_for}.list"])

    return result


rule sss_log_summary:
    input:
        sss_log_input()
    output:
        report(out_dir + "logs/sss_logs.html", caption="../report/Summary_of_SSS_logs.rst")
    conda:
        "../envs/sciL3_env.yaml"
    log:
        stderr = out_dir + f"logs/compile_sss.{aligner}_{assembly}.err",
        stdout = out_dir + f"logs/compile_sss.{aligner}_{assembly}.out"
    shell:
        '''
        python workflow/scripts/compile_sss_logs.py \
        -i {out_dir}logs/ \
        -f 'html' \
        -o {output} > {log.stdout} 2> {log.stderr}
        '''


def bkg_input():
    result = []

    if SSS_FILES == None:
        samples = ALL_SAMPLES
    else:
        samples = SSS_SAMPLES
        
    for s in samples:
        if config['run_filter']:
            result.append([out_dir + f"files_{aligner}_{assembly}/{s}_{summary_for}.list"])
            out_dir + f"files_{aligner}_{assembly}/{{sample}}_filt.list"
        else:
            result.append([out_dir + f"files_{aligner}_{assembly}/{s}_{summary_for}.list"])

    return result


rule bkg_estimate:
    input:
        bkg_input()
    output:
        html = report(out_dir + "logs/background_estimate.html", caption="../report/Background_estimate.rst"),
        csv = report(out_dir + "logs/background_estimate.csv", caption="../report/Background_estimate.rst")
    conda:
        "../envs/r_env.yaml"
    log:
        stderr = out_dir + f"logs/bkg_estimate.err",
        stdout = out_dir + f"logs/bkg_estimate.out"
    threads:
        20
    shell:
        '''
        Rscript workflow/scripts/background_estimate.R \
        -i {out_dir}split_bam_{aligner}_{assembly}/ \
        -w 10000000 \
        -t {threads} \
        -o {output.csv} \
        -o2 {output.html} > {log.stdout} 2> {log.stderr}
        '''


# TODO add Tn5 plate plots