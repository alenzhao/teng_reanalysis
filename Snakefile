SAMPLE_NAMES = ['ENCLB037ZZZ', 'ENCLB038ZZZ', 'ENCLB055ZZZ', 'ENCLB056ZZZ']

TRANSCRIPTOME_FA = 'index/gencode.v16.pc_transcripts.fa'
INDEX = 'index/gencode.v16.pc_transcripts'
INDEX_LIST = expand(INDEX + '.{i}.bt2', i = range(1, 5))

UPDATED_PATH = 'PATH=software/bin:$PATH'

BOWTIE = 'software/bin/bowtie2'
EXPRESS = 'software/bin/express'

BOWTIE_ORIENTATION = ['fr', 'rf', 'ff']
EXPRESS_ORIENTATION = ['fr', 'rf']

rule all:
    input:
        BOWTIE,
        EXPRESS,
        # INDEX_LIST,

        expand('data/{sample}/r{i}.fastq.gz', sample = SAMPLE_NAMES,
            i = [1, 2]),

        expand('alignments/{sample}/{orientation}/bowtie.bam',
            sample = SAMPLE_NAMES, orientation = ['fr', 'rf', 'ff']),

        # use EXPRESS_ORIENTATION because BOWTIE_ORIENTATION (ff)
        # produces an error since there are no reads
        expand('results/{sample}/{or_align}/{or_quant}/express/results.xprs',
            sample = SAMPLE_NAMES, or_align = EXPRESS_ORIENTATION,
            or_quant = EXPRESS_ORIENTATION),

        expand('results/{sample}/{or_align}/express_default/results.xprs',
            sample = SAMPLE_NAMES, or_align = EXPRESS_ORIENTATION),

        expand('results/{sample}/paper/express/results.xprs',
            sample = SAMPLE_NAMES),
        expand('results/{sample}/paper_forward/express/results.xprs',
            sample = SAMPLE_NAMES),
        expand('results/{sample}/paper_reverse/express/results.xprs',
            sample = SAMPLE_NAMES)



rule get_encode_1:
    output:
        'data/ENCLB037ZZZ/r1.fastq.gz',
        'data/ENCLB037ZZZ/r2.fastq.gz'
    params:
        sample = 'ENCLB037ZZZ'
    threads: 1
    shell:
        'wget -O data/{params.sample}/r1.fastq.gz https://www.encodeproject.org/files/ENCFF001REK/@@download/ENCFF001REK.fastq.gz'
        ' && '
        'wget -O data/{params.sample}/r2.fastq.gz https://www.encodeproject.org/files/ENCFF001REJ/@@download/ENCFF001REJ.fastq.gz'

rule get_encode_2:
    output:
        'data/ENCLB038ZZZ/r1.fastq.gz',
        'data/ENCLB038ZZZ/r2.fastq.gz'
    params:
        sample = 'ENCLB038ZZZ'
    threads: 1
    shell:
        'wget -O data/{params.sample}/r1.fastq.gz https://www.encodeproject.org/files/ENCFF001REI/@@download/ENCFF001REI.fastq.gz'
        ' && '
        'wget -O data/{params.sample}/r2.fastq.gz https://www.encodeproject.org/files/ENCFF001REH/@@download/ENCFF001REH.fastq.gz'

rule get_encode_3:
    output:
        expand('data/ENCLB055ZZZ/r{i}.fastq.gz', i = [1, 2])
    params:
        sample = 'ENCLB055ZZZ'
    threads: 1
    shell:
        'wget -O data/{params.sample}/r1.fastq.gz https://www.encodeproject.org/files/ENCFF001RED/@@download/ENCFF001RED.fastq.gz'
        ' && '
        'wget -O data/{params.sample}/r2.fastq.gz https://www.encodeproject.org/files/ENCFF001RDZ/@@download/ENCFF001RDZ.fastq.gz'

rule get_encode_4:
    output:
        expand('data/ENCLB056ZZZ/r{i}.fastq.gz', i = [1, 2])
    params:
        sample = 'ENCLB056ZZZ'
    threads: 1
    shell:
        'wget -O data/{params.sample}/r1.fastq.gz https://www.encodeproject.org/files/ENCFF001REG/@@download/ENCFF001REG.fastq.gz'
        ' && '
        'wget -O data/{params.sample}/r2.fastq.gz https://www.encodeproject.org/files/ENCFF001REF/@@download/ENCFF001REF.fastq.gz'

###
# annotation
###

rule get_transcriptome:
    output:
        TRANSCRIPTOME_FA
    threads: 1
    shell:
        'wget -O {TRANSCRIPTOME_FA}.gz --quiet ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_16/gencode.v16.pc_transcripts.fa.gz'
        ' && '
        'gunzip {TRANSCRIPTOME_FA}.gz'


rule build_index:
    input:
        BOWTIE,
        TRANSCRIPTOME_FA
    output:
        # INDEX + '.1.bt2'
        INDEX_LIST
        # expand(INDEX + '.{i}.bt2', i = range(1, 5)),
        # expand(INDEX + '.rev.{i}.bt2', i = range(1, 2))
    shell:
        '{UPDATED_PATH} bowtie2-build'
        ' --seed 42'
        ' --offrate 1'
        ' {TRANSCRIPTOME_FA}'
        ' {INDEX}'

###
# software
###

rule get_bowtie:
    output:
        BOWTIE
    params:
        url = 'https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip'
    shell:
        'wget -O bowtie.zip {params.url}'
        ' && '
        'unzip bowtie.zip'
        ' && '
        'find bowtie* -perm /111 -type f -exec cp {{}} software/bin \;'
        ' && '
        'rm -rf bowtie*'

rule get_express:
    output:
        EXPRESS
    params:
        url = 'http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz'
    shell:
        'wget -O express.tar.gz {params.url}'
        ' && '
        'tar -xf express.tar.gz'
        ' && '
        'find express* -perm /111 -type f -exec cp {{}} software/bin \;'
        ' && '
        'rm -rf express*'

###
# running the samples
###

rule align_paper:
    input:
        'data/{sample}/r1.fastq.gz',
        'data/{sample}/r2.fastq.gz',
        INDEX_LIST
    output:
        'alignments/{sample}/paper/bowtie.bam'
    benchmark: 'logs/{sample}/bowtie/paper/benchmark.log'
    log: 'logs/{sample}/bowtie/paper/run.log'
    threads: 20
    shell:
        '{UPDATED_PATH} bowtie2'
        ' -a'
        ' -x {INDEX}'
        ' -p {threads}'
        ' -1 {input[0]} '
        ' -2 {input[1]} | '
        'samtools view -Sb - > {output}'

rule align_fr:
    input:
        'data/{sample}/r1.fastq.gz',
        'data/{sample}/r2.fastq.gz',
        INDEX_LIST
    output:
        'alignments/{sample}/fr/bowtie.bam'
    threads: 20
    shell:
        '{UPDATED_PATH} bowtie2'
        ' --fr '
        '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 '
        '--no-mixed --no-discordant '
        '--seed 42 '
        '-k 200 '
        '-x {INDEX} '
        '-p {threads} '
        '-1 {input[0]} '
        '-2 {input[1]} | '
        'samtools view -Sb - > {output}'

rule align_rf:
    input:
        'data/{sample}/r1.fastq.gz',
        'data/{sample}/r2.fastq.gz',
        INDEX_LIST
    output:
        'alignments/{sample}/rf/bowtie.bam'
    threads: 20
    shell:
        '{UPDATED_PATH} bowtie2'
        ' --rf '
        '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 '
        '--no-mixed --no-discordant '
        '--seed 42 '
        '-k 200 '
        '-x {INDEX} '
        '-p {threads} '
        '-1 {input[0]} '
        '-2 {input[1]} | '
        'samtools view -Sb - > {output}'

rule align_ff:
    input:
        'data/{sample}/r1.fastq.gz',
        'data/{sample}/r2.fastq.gz',
        INDEX_LIST
    output:
        'alignments/{sample}/ff/bowtie.bam'
    threads: 20
    shell:
        '{UPDATED_PATH} bowtie2'
        ' --ff '
        '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 '
        '--no-mixed --no-discordant '
        '--seed 42 '
        '-k 200 '
        '-x {INDEX} '
        '-p {threads} '
        '-1 {input[0]} '
        '-2 {input[1]} | '
        'samtools view -Sb - > {output}'

rule express_stranded:
    input:
        'alignments/{sample}/{or_align}/bowtie.bam'
    output:
        'results/{sample}/{or_align}/{or_quant}/express/results.xprs'
    benchmark: 'logs/{sample}/express/{or_align}/{or_quant}/benchmark.log'
    log: 'logs/{sample}/express/{or_align}/{or_quant}/run.log'
    threads: 2
    params:
        out = 'results/{sample}/{or_align}/{or_quant}/express'
    shell:
        '{UPDATED_PATH} express'
        ' --{wildcards.or_quant}-stranded'
        ' -o {params.out}'
        ' {TRANSCRIPTOME_FA}'
        ' {input}'

rule express_default:
    input:
        'alignments/{sample}/{or_align}/bowtie.bam'
    output:
        'results/{sample}/{or_align}/express_default/results.xprs'
    benchmark: 'logs/{sample}/express_default/{or_align}/benchmark.log'
    log: 'logs/{sample}/express_default/{or_align}/run.log'
    threads: 2
    params:
        out = 'results/{sample}/{or_align}/express_default'
    shell:
        '{UPDATED_PATH} express'
        ' -o {params.out}'
        ' {TRANSCRIPTOME_FA}'
        ' {input}'

rule express_paper:
    input:
        'alignments/{sample}/paper/bowtie.bam'
    output:
        'results/{sample}/paper/express/results.xprs'
    benchmark: 'logs/{sample}/express/paper/benchmark.log'
    log: 'logs/{sample}/express/paper/run.log'
    threads: 2
    params:
        out = 'results/{sample}/paper/express'
    shell:
        '{UPDATED_PATH} express'
        ' -o {params.out}'
        ' {TRANSCRIPTOME_FA}'
        ' {input}'

rule express_paper_reverse:
    input:
        'alignments/{sample}/paper/bowtie.bam'
    output:
        'results/{sample}/paper_reverse/express/results.xprs'
    benchmark: 'logs/{sample}/express/paper_reverse/benchmark.log'
    log: 'logs/{sample}/express/paper_reverse/run.log'
    threads: 2
    params:
        out = 'results/{sample}/paper_reverse/express'
    shell:
        '{UPDATED_PATH} express'
        ' --rf-stranded'
        ' -o {params.out}'
        ' {TRANSCRIPTOME_FA}'
        ' {input}'

rule express_paper_forward:
    input:
        'alignments/{sample}/paper/bowtie.bam'
    output:
        'results/{sample}/paper_forward/express/results.xprs'
    benchmark: 'logs/{sample}/express/paper_forward/benchmark.log'
    log: 'logs/{sample}/express/paper_forward/run.log'
    threads: 2
    params:
        out = 'results/{sample}/paper_forward/express'
    shell:
        '{UPDATED_PATH} express'
        ' --fr-stranded'
        ' -o {params.out}'
        ' {TRANSCRIPTOME_FA}'
        ' {input}'
