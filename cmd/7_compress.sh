#!/usr/bin/env bash
set -euo pipefail

echo "----- 7_compress.sh -----"
echo "Merging ..."
## ----------------------------------------------------------------------------------------------------
find nucl/clustered -maxdepth 1 -name '*rep_seq.fasta' -exec cat {} \; > nucl/nucl_c.fa
cat nucl/nucl_*.fa | seqkit sort --quiet | seqkit shuffle --quiet -s 0 > nucl/nucl.fa

python -c "
import pandas as pd
from collections import defaultdict

sequence = defaultdict(list)
accession = set()
with open('nucl/nucl.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line.split()
            accession.add(ls[0][1:].rsplit('_', 1)[0])
            subset = ls[1].split('|')[1]
        sequence[subset.replace('/', ';')].append(line)

for key, val in sequence.items():
    with open('nucl/sarg.{}.fa'.format(key), 'w') as w:
        w.write(''.join(val))

metadata = pd.read_table('nucl/raw.tsv')
metadata[metadata.accession.isin(accession)].sort_values('accession').to_csv('nucl/sarg.metadata.tsv', index=False, sep='\t')
"
## ----------------------------------------------------------------------------------------------------





echo "Compressing ..."
## ----------------------------------------------------------------------------------------------------
wget -qN https://zenodo.org/records/19924976/files/database.tar.gz
tar -zxf database.tar.gz

cp prot/prot.fa database/sarg.fa
cp nucl/sarg.metadata.tsv nucl/sarg*.fa database

tar --sort=name -zcf database.tar.gz database
## ----------------------------------------------------------------------------------------------------
