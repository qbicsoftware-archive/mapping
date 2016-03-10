import os
import sys
import subprocess
import tempfile
import uuid
import shutil
from datetime import datetime
from os.path import join as pjoin
from os.path import exists as pexists
from xml.etree import ElementTree
import hashlib
import base64
import csv
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


rule merge_groups:
    input: expand("map_bwa/{group}.bam", group=GROUPS['group'])
    output: result("aligned.bam")
    shell: "samtools merge -l 9 -@ 5 {output} {input}"


rule merge_sam_index:
    input:  "{name}.bam"
    output: "{name}.bam.bai"
    shell: "samtools index {input}"
