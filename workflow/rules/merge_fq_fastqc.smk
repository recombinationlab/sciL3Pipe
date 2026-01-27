

rule merge_fastqs_pe:
    '''
    Rename and if needed merge fastq files
    '''
    input: 
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: 
        r1 = out_dir + "merged/{sample}.R1.fastq.gz",
        r2 = out_dir + "merged/{sample}.R2.fastq.gz"
    conda:
        "../envs/sciL3_env.yaml"
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


rule fastqc:
    input:
        [out_dir + "merged/{sample}.R1.fastq.gz",
         out_dir + "merged/{sample}.R2.fastq.gz"]
        # [lambda wildcards: FILES[wildcards.sample]['R1'],
        # lambda wildcards: FILES[wildcards.sample]['R2']]
    output:
        [out_dir + "fastqc/{sample}.R1_fastqc.html",
         out_dir + "fastqc/{sample}.R2_fastqc.html"]
    conda:
        "../envs/sciL3_env.yaml"
    params:
        fqc_out_dir = out_dir + "fastqc/"
    log:
        stdout = out_dir + "logs/{sample}_fastqc.out",
        stderr = out_dir + "logs/{sample}_fastqc.err"
    shell:
        '''
        fastqc -o {params.fqc_out_dir} {input} > {log.stdout} 2> {log.stderr}
        '''