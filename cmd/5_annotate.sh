#!/usr/bin/env bash
set -euo pipefail

echo "----- 5_annotate.sh -----"
echo "Retrieving plasmid sequences ..."
## ----------------------------------------------------------------------------------------------------
python -c "
import glob
import gzip
import multiprocessing
multiprocessing.set_start_method('fork', force=True)

from Bio import SeqIO
from tqdm.contrib.concurrent import process_map

def parser(file):
    records = []
    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if ' plasmid' in record.description or ' megaplasmid' in record.description:
                records.append(record)
    return records

r = process_map(parser, glob.glob('plasmid/fna/*/*.fna.gz'), max_workers=32, chunksize=1, leave=False)
with open('plasmid/plasmid.fna', 'w') as output_handle:
    SeqIO.write([x for y in r for x in y], output_handle, 'fasta')
"
## ----------------------------------------------------------------------------------------------------





echo "Annotating ARGs on plasmid sequences ..."
## ----------------------------------------------------------------------------------------------------
diamond blastx \
    --db prot/prot.fa \
    --query plasmid/plasmid.fna \
    --out plasmid/plasmid_sarg.txt \
    --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
    --evalue 1e-25 --subject-cover 90 --id 90 \
    --range-culling --frameshift 15 --range-cover 25 \
    --max-hsps 0 --max-target-seqs 25 \
    --threads 32 --quiet

python -c "
import xml.etree.ElementTree as ET

from collections import defaultdict
from math import floor, ceil
from tqdm import tqdm
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

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

ids = set()
with open('plasmid/plasmid_sarg.txt') as f:
    for line in f:
        ls = line.rstrip().split()
        qlen, qstart, qend = int(ls[4]), int(ls[5]), int(ls[6])
        if max(qstart, qend) + 5000 > qlen or min(qstart, qend) - 5000 < 0:
            ids.add(ls[0])

qseqid2topology = {x[0]: x[1] for x in get_topology(ids)}

bed = defaultdict(list)
qrange = defaultdict(set)
with open('plasmid/plasmid_sarg.txt') as f:
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
with open('plasmid/plasmid.fna') as handle:
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
                records.append(SeqRecord(subseq, id=f'plasmid@{record.id}_{row[0]+1}-{row[1]}:{row[2]}', description = row[3]))

with open('plasmid/plasmid_arg.fa', 'w') as output_handle:
    SeqIO.write(records, output_handle, 'fasta')
"
## ----------------------------------------------------------------------------------------------------





echo "Filtering plasmid sequences ..."
## ----------------------------------------------------------------------------------------------------
seqkit grep --quiet -f <(cut -f1 plasmid/plasmid_sarg.txt) plasmid/plasmid.fna > plasmid/plasmid_sub.fa
genomad end-to-end --cleanup plasmid/plasmid_sub.fa plasmid/genomad genomad_db \
    --disable-find-proviruses \
    --conservative \
    --threads 32 --quiet

genomad end-to-end --cleanup plasmid/plasmid_arg.fa plasmid/genomad genomad_db \
    --disable-find-proviruses \
    --relaxed \
    --threads 32 --quiet

python -c "
import pandas as pd
from Bio import SeqIO

fea = pd.read_table('plasmid/genomad/plasmid_sub_marker_classification/plasmid_sub_features.tsv')
sub = pd.read_table('plasmid/genomad/plasmid_sub_summary/plasmid_sub_plasmid_summary.tsv')
sub = set(sub[sub.seq_name.isin(fea[fea.n_uscg == 0].seq_name)].seq_name)

arg = pd.read_table('plasmid/genomad/plasmid_arg_aggregated_classification/plasmid_arg_aggregated_classification.tsv')
arg = arg[arg.seq_name.str.split('@').str.get(-1).str.rsplit('_', n=1).str.get(0).isin(sub)]
arg = arg[arg.plasmid_score / arg.chromosome_score > 2]
arg = set(arg.seq_name)

records = []
with open('plasmid/plasmid_arg.fa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        if record.id in arg:
            records.append(record)

with open('nucl/seq/plasmid.fa', 'w') as output_handle:
    SeqIO.write(records, output_handle, 'fasta')
"
## ----------------------------------------------------------------------------------------------------
