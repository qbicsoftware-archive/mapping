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
INI_PATH = config['etc']


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

OUTPUT_FILES = []

OUTPUT_FILES.extend(expand(result("FastQC_{name}.zip"), name=INPUT_FILES))
OUTPUT_FILES.extend(expand(result("FastQC_{name}.html"), name=INPUT_FILES))

rule all:
    input: OUTPUT_FILES

rule checksums:
    output: "checksums.ok"
    run:
        out = os.path.abspath(str(output))
        shell("cd %s; "
              "sha256sum -c *.sha256sum && "
              "touch %s" % (data('.'), out))

rule flagstat:
    input: "{file}.sorted.bam"
    output: "{file}.sorted.bam.flagstat"
    shell: "samtools flagstat {input} > {output}"

rule fastqc:
    input: data("{name}.fastq.gz")
    output: "fastqc/{name}"
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        #TODO jobscript load module
        shell("fastqc {input} -o {output}")


def files_by_lane(wildcards):
  lane = wildcards["lane"]
  group = wildcards["group"]
  return pares()


rule map_bwa:
    input: files_by_lane
    output: "map_bwa/{group}_{lane}.bam"
    threads: 10
    run:
        # TODO additional header LB (library)?
        options = r"-M -R '@RG\tID:{group}\tSM:{group}'".format(**wildcards)
        options = options + " " + ref(str(config['params']['fasta']))
        shell("bwa mem -t {threads} %s {input}/forward.fastq.gz {input}/backward.fastq.gz | samtools view -Sb - | samtools sort -O bam -T $(mktemp) - > {output}" % options)

rule sam_sort:
    input: "{file}.unsorted.bam"
    output: "{file}.sorted.bam"
    shell: "samtools sort -O bam -T $(mktemp) -o {output} {input}"

rule sam_index:
    input: "{file}.sorted.bam"
    output: "{file}.sorted.bam.bai"
    run:
        shell("samtools index {input}")

rule merge_groups:
    input: "map_bwa/{pool}_{group}.bam"
    output: result("merged/aligned.bam")
    shell: "samtools merge -l 9 -@ 5 {output} {input}"

rule merge_sam_index:
    input:  result("merged_bam/aligned.bam")
    output: result("merged_bam/aligned.bam.bai")
    run:
        shell("samtools index {input}")