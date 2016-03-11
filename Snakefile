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


GROUPS = pd.read_csv(etc('groups.txt'), sep='\t')


OUTPUT_FILES = []
OUTPUT_FILES.extend(expand(result("fastqc/{name}"), name=INPUT_FILES))

rule all:
    input: result("aligned.bam.bai"), OUTPUT_FILES, "checksums.ok"

rule checksums:
    output: "checksums.ok"
    run:
        out = os.path.abspath(str(output))
        shell("cd %s; "
              "sha256sum -c *.sha256sum && "
              "touch %s" % (data('.'), out))

rule fastqc:
    input: data("{name}.fastq.gz")
    output: result("fastqc/{name}")
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        #TODO jobscript load module
        shell("fastqc {input} -o {output}")


def files_by_group(wildcards):
    group = wildcards["group"]
    files = GROUPS[GROUPS['group'] == group]
    assert len(files) > 0
    first = files[files['file'].str.contains('_R1')]
    first = first.sort_values(by='file').reset_index(drop=True)
    second = files[files['file'].str.contains('_R2')]
    second = second.sort_values(by='file').reset_index(drop=True)
    assert(len(files) == len(first) + len(second))
    if not first['file'].str.replace('_R1', '_R2').equals(second['file']):
        raise ValueError("Invalid file names for paired end sequencing")
    first = list(sorted(first['file'].values))
    second = list(sorted(second['file'].values))
    return [data(file) for file in first + second]


rule bwa_mem:
    input: files_by_group
    output: "map_bwa/{group}.bam"
    params: r"-t 10 -M -R '@RG\tID:{group}\tSM:{group}'"
    run:
        fasta = ref(config['params']['fasta'])
        first = " ".join(input[:len(input) / 2])
        second = " ".join(input[len(input) / 2:])
        assert len(first) == len(second)
        shell(
            "bwa mem {params} %s "
                "<(cat %s) <(cat %s) |" % (fasta, first, second) +
            "samtools sort -o -l1 -@2 -T $(mktemp) - > {output}"
        )


def bamfile_by_group(wildcards):
    group = wildcards['group']
    file = GROUPS[GROUPS['group'] == group]
    assert len(file) == 1
    return [data(file['file'].values[0]), data(file['file'].values[0] + '.bai')]


rule bwa_mem_bam:
    input: bamfile_by_group
    output: "map_bwa/{group}.bam"
    params: "-t10", "-M", "-R", r"@RG\tID:{group}\tSM:{group}"
    run:
        fasta = ref(config['params']['fasta'])
        with tempfile.TemporaryDirectory() as tmp:
            sorted_bam = os.path.join(tmp, 'sorted.bam')
            first = os.path.join(tmp, 'first.fastq')
            mapped = os.path.join(tmp, 'mapped')
            sort_temp = os.path.join(tmp, 'sort_tmp')
            sort_temp2 = os.path.join(tmp, 'sort_tmp2')
            os.mkdir(sort_temp)
            os.mkfifo(sorted_bam)
            os.mkfifo(first)
            os.mkfifo(mapped)
            sort = subprocess.Popen(['samtools', 'sort', '-n', '-@', '5',
                                     '-o', sorted_bam, '-O', 'bam',
                                     '-T', sort_temp2,
                                     str(input[0])])
            to_fastq = subprocess.Popen(
                ['bedtools', 'bamtofastq', '-i', sorted_bam, '-fq', first]
            )
            mapping = subprocess.Popen(
                ['bwa', 'mem', '-p'] + list(params) + [fasta, first],
                stdout=subprocess.PIPE,
            )
            sort = subprocess.Popen(
                ["samtools", "sort", "-l", "1", "-@", "5", "-T", sort_temp,
                 '-O', 'bam', "-o", str(output), '-'],
                stdin=mapping.stdout,
            )
            mapping.stdout.close()
            retcode = sort.wait()
            assert retcode == 0


rule merge_groups:
    input: expand("map_bwa/{group}.bam", group=GROUPS['group'])
    output: result("aligned.bam")
    shell: "samtools merge -l 9 -@ 5 {output} {input}"


rule bam_index:
    input:  "{name}.bam"
    output: "{name}.bam.bai"
    shell: "samtools index {input}"
