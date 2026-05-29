#!/usr/bin/env bash
set -euo pipefail

echo "----- 4_annotate.sh -----"
echo "Annotating ARGs ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/out

find assembly/fna -maxdepth 1 -name '*.fna.gz' | sort | xargs -P 4 -I {} bash -c '
    filename=${1%.fna*};
    filename=${filename##*/};
    diamond blastx \
        --db prot/prot.fa \
        --query $1 \
        --out nucl/out/${filename}.txt \
        --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
        --evalue 1e-25 --subject-cover 90 --id 90 \
        --range-culling --frameshift 15 --range-cover 25 \
        --max-hsps 0 --max-target-seqs 25 \
        --threads 8 --quiet
' - {}
## ----------------------------------------------------------------------------------------------------





echo "Extracting sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/seq 

python -c "
import pandas as pd
import glob
import gzip
import os
import xml.etree.ElementTree as ET
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

from collections import defaultdict
from math import floor, ceil
from tqdm import tqdm
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm.contrib.concurrent import process_map

Entrez.email = 'A.N.Other@example.com'

def sort_coordinate(start, end):
    return (start - 1, end, '+') if start < end else (end - 1, start, '-')

def compute_overlap(coordinates):
    qstart, qend, sstart, send = coordinates
    overlap = min(qend, send) - max(qstart, sstart)
    return max(overlap / (qend - qstart), overlap / (send - sstart))

def get_topology(ids, chunksize=1000):
    topology = []
    n = len(ids) // chunksize + 1
    folds = [list(ids)[i::n] for i in range(n)]
    for fold in tqdm(folds, leave=False):
        tmp = []
        while len(tmp) != len(fold):
            with Entrez.esummary(db='nuccore', id=','.join(fold), rettype='xml', version='2.0') as handle:
                root = ET.parse(handle).getroot()
                tmp = [[i,j] for i,j in zip([x.text for x in root.findall('.//AccessionVersion')], [x.text for x in root.findall('.//Topology')])]

            try: assert len(tmp) == len(fold)
            except: 
                todo = set(fold) - set(x[0] for x in tmp)
                with Entrez.efetch(db='nuccore', id=','.join(todo), rettype='gb', retmode='xml') as handle:
                    for i in Entrez.read(handle):
                        tmp.append([i['GBSeq_accession-version'], i['GBSeq_topology']])
        topology.extend(tmp)

    return topology

def extract_sequence(file):
    filename = os.path.basename(file).split('.fna.gz')[0]

    bed = defaultdict(list)
    qrange = defaultdict(set)
    with open(f'nucl/out/{filename}.txt') as f:
        for line in f:
            ls = line.split()
            qseqid, sseqid, qlen = ls[0], ls[1], int(ls[4])
            qstart, qend, strand = sort_coordinate(int(ls[5]), int(ls[6]))
            if (
                qseqid not in qrange or
                all([compute_overlap((qstart, qend, *x)) < 0.25 for x in qrange.get(qseqid)])
            ):
                qrange[qseqid].add((qstart, qend))
                qcoord = floor((qstart + qend) / 2) if strand == '+' else ceil((qstart + qend) / 2)

                topology = qseqid2topology.get(qseqid)
                if (qend + 2500 < qlen and qstart - 2500 > 0) or topology == 'circular':
                    bed[ls[0]].append((qcoord - 5000, qcoord + 5000, strand, sseqid, topology))

    records = []
    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in bed:
                for row in bed.get(record.id):
                    subseq = record.seq[max(row[0], 0): min(row[1], len(record.seq))]
                    fullseq = record.seq * (5000 // len(record.seq) + 1)

                    if row[-1] == 'circular':
                        if row[0] < 0:
                            subseq = fullseq[(len(record.seq) + row[0]):] + subseq
                        if row[1] > len(record.seq):
                            subseq = subseq + fullseq[:(row[1] - len(record.seq))]

                    if row[2] == '-':
                        subseq = subseq.reverse_complement()
                    records.append(SeqRecord(subseq, id=f'{record.id}_{row[0]+1}-{row[1]}:{row[2]}', description = row[3]))

    with open(f'nucl/seq/{filename}.fa', 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

gtdb = pd.concat([
    pd.read_table('assembly/bac120_metadata.tsv'),
    pd.read_table('assembly/ar53_metadata.tsv')
])
accession = set(gtdb[gtdb['ncbi_assembly_level'].isin({'Chromosome', 'Complete Genome'})].accession.str.split('_', n=1).str.get(-1))

aset = set()
with open('assembly/accession2assembly.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        if ls[1] in accession:
            aset.add(ls[0])

ids = set()
for file in tqdm(glob.glob('nucl/out/*.txt'), leave=False):
    with open(file) as f:
        for line in f:
            ls = line.rstrip().split()
            qlen, qstart, qend = int(ls[4]), int(ls[5]), int(ls[6])
            if max(qstart, qend) + 5000 > qlen or min(qstart, qend) - 5000 < 0:
                if ls[0] in aset:
                    ids.add(ls[0])

qseqid2topology = {x[0]: x[1] for x in get_topology(ids)}
process_map(extract_sequence, glob.glob('assembly/fna/*.fna.gz'), max_workers=32, chunksize=1, leave=False)
"
## ----------------------------------------------------------------------------------------------------
