#!/usr/bin/env bash
set -euo pipefail

echo "----- 3_download.sh -----"
echo "Downloading RefSeq plasmids ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p plasmid/fna

python -c "
import pandas as pd

assembly = pd.read_table('assembly/assembly_summary_refseq.txt', skiprows=1, low_memory=False).rename({'#assembly_accession': 'assembly'}, axis=1)
assembly = assembly[(assembly.group.isin(['archaea', 'bacteria'])) & ((assembly['ftp_path'] != 'na'))]
assembly = assembly[assembly.assembly_level.isin({'Complete Genome', 'Chromosome'})]
assembly['fna'] = assembly['ftp_path'] + assembly['ftp_path'].str.split('/').str.get(-2) + '_genomic.fna.gz'

for kingdom in ['archaea', 'bacteria']:
    fna = assembly[assembly['group'] == kingdom]['fna'].to_list()
    n = 2 if kingdom == 'archaea' else 32
    chunks = [fna[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        with open('plasmid/fna/' + kingdom + '.split.' + str(i) + '.id', 'w') as w:
            w.write('\n'.join(chunk) + '\n')

with open('plasmid/assertion.tsv', 'w') as f:
    f.write(f'{assembly.assembly.nunique()}\n')
"

find plasmid/fna -maxdepth 1 -name '*.id' | sort | xargs -P 8 -I {} bash -c '
    wget -i ${1} -qN --retry-on-http-error=503 -P ${1%.id};
    find ${1%.id} -maxdepth 1 -name "*.fna.gz" | sort | xargs cat > ${1%.id}.fna.gz;
' - {}

expected_count=$(head -n 1 plasmid/assertion.tsv)
file_count=$(find plasmid/fna -mindepth 2 -name '*.fna.gz' | wc -l | xargs)
if (( file_count != expected_count )); then
    echo "Assertion failed: Found $file_count files, but expected $expected_count. An error may have occurred; please check whether any assemblies are deprecated or rerun this script."
fi
## ----------------------------------------------------------------------------------------------------





echo "Downloading SARG+ ..."
## ----------------------------------------------------------------------------------------------------
mkdir -p prot
wget -q -O prot/prot.fa https://zenodo.org/records/20373707/files/sarg.fa
## ----------------------------------------------------------------------------------------------------
