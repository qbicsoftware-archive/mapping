import os
import sys
import subprocess
import tempfile
from datetime import datetime
from os.path import join as pjoin
from os.path import exists as pexists
import glob
import numpy as np
import pandas as pd
import time

configfile: "config.json"
workdir: config['var']

SNAKEDIR = config['src']

try:
    VERSION = subprocess.check_output(
        ['git', 'describe', '--tags', '--always', '--dirty'],
        cwd=SNAKEDIR
    ).decode().strip()
except subprocess.CalledProcessError:
    VERSION = 'unknown'

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
ETC = config['etc']


def data(path):
    return os.path.join(DATA, path)

def ref(path):
    return os.path.join(REF, path)

def log(path):
    return os.path.join(LOGS, path)

def result(path):
    return os.path.join(RESULT, path)

def etc(path):
    return os.path.join(ETC, path)

if 'params' not in config:
    config['params'] = {}

INPUT_FILES = []
for name in os.listdir(DATA):
    if name.lower().endswith('.sha256sum'):
        continue
    if name.lower().endswith('.fastq'):
        if not name.endswith('.fastq'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-6])
    elif name.lower().endswith('.fastq.gz'):
        if not name.endswith('.fastq.gz'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-len('.fastq.gz')])
    elif name.lower().endswith('.bam'):
        if not name.endswith('.bam'):
            print("Extension bam is case sensitive.", file=sys.stderr)
            exit(1)
    elif name.endswith('.bam.bai'):
        continue
    else:
        print("Unknown data file: %s" % name, file=sys.stderr)
        exit(1)

if len(set(INPUT_FILES)) != len(INPUT_FILES):
    print("Some input file names are not unique")
    exit(1)


DESIGN = pd.read_csv(etc('GROUPS'), sep='\t')

RUNS = {}
for group, df in DESIGN.groupby('group'):
    first = df[df['file'].str.contains('_R1')]
    first = first.sort_values(by='file').reset_index(drop=True)
    second = df[df['file'].str.contains('_R2')]
    second = second.sort_values(by='file').reset_index(drop=True)
    assert(len(df) == len(first) + len(second))
    if not first['file'].str.replace('_R1', '_R2').equals(second['file']):
        raise ValueError("Invalid file names for paired end sequencing")
    first = list(sorted(first['file'].values))
    second = list(sorted(second['file'].values))
    prefix = [name[:name.find('_R1')] for name in first]
    for prefix, file1, file2 in zip(prefix, first, second):
        assert (group, prefix) not in RUNS
        assert prefix and group
        assert '___' not in group and '___' not in prefix
        assert file1.startswith(prefix) and file2.startswith(prefix)
        RUNS[(group, prefix)] = (file1, file2)



OUTPUT_FILES = []
OUTPUT_FILES.extend(expand(result("fastqc/{name}"), name=INPUT_FILES))

rule all:
    input: result("aligned.bam.bai"), OUTPUT_FILES, "checksums.ok"

rule checksums:
    output: "checksums.ok"
    threads: 1
    run:
        out = os.path.abspath(str(output))
        cksums = glob.glob(data("*.sha256sum"))
        if cksums:
            shell("cd %s; "
                  "sha256sum -c *.sha256sum && "
                  "touch %s" % (data('.'), out))
        else:
            shell("touch %s" % out)

rule fastqc:
    input: data("{name}.fastq")
    output: result("fastqc/{name}")
    threads: 1
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        shell("fastqc {input} -o {output}")


def align_sort(fastq, outfile, fasta, tmp, params=[]):
    if not isinstance(fastq, (tuple, list)):
        fastq = [fastq]
    assert len(fastq) in [1, 2]

    if len(fastq) == 1:
        params = ['-p'] + params
    bwa = subprocess.Popen(
        ['bwa', 'mem'] + params + [fasta] + list(fastq),
        stdout=subprocess.PIPE,
        stderr=sys.stderr
    )

    sort_tmp = os.path.join(tmp,'sortedtmp.bam')

    sort = subprocess.Popen(
        ['samtools', 'sort', '-l', '0', '-@', '10',
         '-T', os.path.join(tmp, 'sort_tmp'),
         '-O', 'bam', '-o', sort_tmp, '-'],
        stdin=bwa.stdout,
        stdout=sys.stdout,
        stderr=sys.stderr
    )

    bwa.stdout.close()
    
    retcodebwa = bwa.wait()
    retcodesort = sort.wait()
    assert retcodebwa == 0
    assert retcodesort == 0
    assert bwa.returncode == 0
    assert sort.returncode == 0
    
    rmdup = subprocess.Popen(['picard','MarkDuplicates','I='+sort_tmp,'O='+str(outfile), 'M='+os.path.join(tmp, 'output_matrix'),'AS=true','VALIDATION_STRINGENCY=LENIENT'],
        stderr=sys.stderr,
        stdout=sys.stdout
    )  

    retcodermdup = rmdup.wait()
    assert retcodermdup == 0
    assert rmdup.returncode ==0

def pipe_bam_to_fastq(bam, fastq, temp_prefix):
    collate = subprocess.Popen(
        ['samtools', 'collate', '-u', '-O', '-T', temp_prefix, bam],
        stdout=subprocess.PIPE,
        stderr=sys.stderr
    )
    to_fastq = subprocess.Popen(
        ['samtools', 'fastq', '-i', '-', fastq],
        stdin=collate.stdout,
        stderr=sys.stderr
    )
    collate.stdout.close()
    return collate, to_fastq

rule trim:
    input: lambda w: [data(p) for p in RUNS[(w['group'], w['prefix'])]]
    output:
       L="trim/{group}___{prefix}_R1.fastq",
       R="trim/{group}___{prefix}_R2.fastq"
    params: "-match_perc 80", "-qcut 15"
    run:
        # TODO more than two elements per group
        # TODO get adapters from fastcq
        # TODO Default adapter is used
        if len(input) == 1:
                assert input[0].endswith('.bam') or input[0].endswith('.cram')
                params.append('-p')
                fastq = [os.path.join(tmp, 'as_fastq.fastq')]
                tmp_prefix = os.path.join(tmp, 'collate')
                progs.extend(pipe_bam_to_fastq(input[0], fastq[0], tmp_prefix))
        elif len(input) == 2:
                assert all((name.endswith('.fastq') or
                            name.endswith('.fastq.gz') or
                            name.endswith('.fq') or
                            name.endswith('.fq.gz')) for name in input)
                fastq = list(input)
        else:
                raise ValueError('Invalid input: %s' % input)
        print(fastq[0])
        print(fastq[1])
        shell("SeqPurge -in1 fastq[0] -in2 fastq[1] "
              + "-out1 output['L'] -out2 output['R'] "
              + params)

rule bwa_mem:
    #input: lambda w: ["trim/" + p for p in RUNS[(w['group'], w['prefix'])]]
    input: L="trim/{group}___{prefix}_R1.fastq",
           R="trim/{group}___{prefix}_R2.fastq"
    output: "map_bwa/{group}___{prefix}.bam"
    params: "-t 20", "-M", "-R", r"@RG\tID:{group}\tSM:{group}"
    run:
        print(input)
        fasta = ref(config['params']['fasta'])
        with tempfile.TemporaryDirectory() as tmp:
            progs = []
            # We accept either 1 bam or 2 fastq files
            if len(input) == 1:
                assert input[0].endswith('.bam') or input[0].endswith('.cram')
                params.append('-p')
                fastq = [os.path.join(tmp, 'as_fastq.fastq')]
                tmp_prefix = os.path.join(tmp, 'collate')
                progs.extend(pipe_bam_to_fastq(input[0], fastq[0], tmp_prefix))
            elif len(input) == 2:
                assert all((name.endswith('.fastq') or
                            name.endswith('.fastq.gz') or
                            name.endswith('.fq') or
                            name.endswith('.fq.gz')) for name in input)
                fastq = list(input)
                print(fastq)
            else:
                raise ValueError('Invalid input: %s' % input)
            align_sort(fastq, output, fasta, tmp, params)
            time.sleep(1)
            for prog in progs:
                assert prog.returncode == 0

rule merge_groups:
    input:
        ["map_bwa/{group}___{prefix}.bam".format(group=key[0], prefix=key[1])
         for key in RUNS]
    output: result("aligned.bam")
    threads: 10
    shell: "samtools merge -l 9 -@ 20 {output} {input}"


rule bam_index:
    input:  "{name}.bam"
    output: "{name}.bam.bai"
    threads: 1
    shell: "samtools index {input}"
