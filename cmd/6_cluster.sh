#!/usr/bin/env bash
set -euo pipefail

echo "----- 6_cluster.sh -----"
echo "Splitting ARG-containing sequences ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/raw
find nucl/seq -maxdepth 1 -name '*.fa' | sort | xargs cat > nucl/raw.fa

python -c "
import pandas as pd
from collections import defaultdict

assembly2species = {}
with open('assembly/assembly2species.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        assembly2species[ls[0]] = ls[1:]

accession2assembly = {}
with open('assembly/accession2assembly.tsv') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        accession2assembly[ls[0]] = ls[1]

accession = set()
sequence = defaultdict(list)
with open('nucl/raw.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line[1:].split()
            if 'plasmid@' in ls[0]:
                lineage = 'plasmid'
            else:
                accession.add(ls[0].rsplit('_', 1)[0])
                lineage = str(assembly2species.get(accession2assembly.get(ls[0].rsplit('_', 1)[0]))[-1])
            subset = ls[1].split('|')[1] + '@' + ls[1].split('|')[2]
        sequence[subset.replace('/','_SLASH_').replace('*','_STAR').replace('''\'''', '_PRIME')+ '.' + lineage].append(line)

with open('nucl/nucl_a.fa', 'w') as v, open('nucl/nucl_b.fa', 'w') as w:
    for key, val in sequence.items():
        record = ''.join(val)
        if record.count('>') == 1:
            if 'plasmid@' in record:
                v.write(record)
            else:
                w.write(record)
        else:
            with open('nucl/raw/{}.fa'.format(key), 'w') as u:
                u.write(record)

with open('nucl/raw.tsv', 'w') as f:
    f.write('\t'.join(['accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) + '\n')
    for i in accession:
        f.write('\t'.join([i] + assembly2species.get(accession2assembly.get(i))[0].split(';')) + '\n')
"
## ----------------------------------------------------------------------------------------------------





echo "Clustering ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p nucl/clustered

find nucl/raw -maxdepth 1 -name '*.fa' | xargs -P 8 -I {} bash -c '
    filename=${1%.fa*};
    filename=${filename##*/};
    mkdir -p nucl/clustered/$filename;
    mmseqs easy-cluster \
        $1 nucl/clustered/$filename/$filename nucl/clustered/$filename/$filename \
        -c 0.9995 --min-seq-id 0.9995 --cov-mode 1 \
        -s 7.5 --cluster-reassign --threads 4 -v 0 > /dev/null;
    mv nucl/clustered/$filename/${filename}_rep_seq.fasta nucl/clustered/${filename}_rep_seq.fasta;
    rm -rf nucl/clustered/$filename/;
' - {}

python -c "
import glob

a = {x.split('/')[-1].split('_rep_seq.fasta')[0] for x in glob.glob('nucl/clustered/*_rep_seq.fasta')}
b = {x.split('/')[-1].split('.fa')[0] for x in glob.glob('nucl/raw/*.fa')}

assert a==b, 'Some files are not processed for some reasons. Try again.'
"
## ----------------------------------------------------------------------------------------------------
