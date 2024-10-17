# Argo-supplementary
Instruction on building the database of Argo.

## Prerequisite
### Step 1: Install necessary packages
```bash
conda install -c bioconda -c conda-forge 'seqkit>=2.6.1' 'genomad>=1.8.1' 'mmseqs2>=15.6f452' 'diamond==2.1.8' 'tqdm' 'pandas' 'biopython'
genomad download-database .
```

### Step 2: Collect GTDB assemblies
> [!NOTE]
> Some genomes (including representative genomes) may be deprecated by NCBI at the time of GTDB release. The representative ones can be downloaded from the FTP of GTDB (`genomic_files_reps/gtdb_genomes_reps.tar.gz`) and manually included. The other non-representative ones cannot be recovered.

```bash
## download metadata files
mkdir -p assembly/fna
wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -P assembly
wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt -P assembly
wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -P assembly
wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt -P assembly

wget -qN --show-progress https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz -P assembly
wget -qN --show-progress https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz -P assembly

gzip -d assembly/ar53_metadata.tsv.gz
gzip -d assembly/bac120_metadata.tsv.gz

## read metadata, split assemblies into chunks according to taxonomy
python -c "
import pandas as pd

## merge gtdb with ncbi ftp link
ncbi = pd.concat([
    pd.read_table('assembly/assembly_summary_refseq.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly/assembly_summary_refseq_historical.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly/assembly_summary_genbank.txt', skiprows=1, low_memory=False),
    pd.read_table('assembly/assembly_summary_genbank_historical.txt', skiprows=1, low_memory=False),
]).rename({'#assembly_accession': 'assembly'}, axis=1)

gtdb = pd.concat([
    pd.read_table('assembly/ar53_metadata.tsv', low_memory=False),
    pd.read_table('assembly/bac120_metadata.tsv', low_memory=False)
])
gtdb['assembly'] = gtdb.accession.str.split('_', n=1).str.get(-1)

## some may no longer be available
assembly = pd.merge(gtdb, ncbi[['assembly', 'ftp_path']], how='left', on='assembly')
assembly[(assembly.ftp_path == 'na') | (assembly.ftp_path.isnull())].to_csv('assembly/deprecated.tsv', index=False, sep='\t')

assembly = assembly[(~(assembly.ftp_path == 'na') | (assembly.ftp_path.isnull())) | (assembly.gtdb_representative == 't')]
assembly['taxonomy'] = assembly['gtdb_taxonomy'].str.replace('[a-z]__', '', regex=True)
assembly['group'] = assembly['taxonomy'].str.split(';').str.get(0).str.lower()

## create an index file for tracing
assembly['index'] = pd.Categorical(assembly['taxonomy'], ordered=False).codes + 1
assembly[['assembly', 'taxonomy', 'index']].to_csv('assembly/assembly2species.tsv', sep='\t', index=False, header=None)

## create a list for wget
assembly['fna'] = assembly['ftp_path'] + '/' + assembly['ftp_path'].str.split('/').str.get(-1) + '_genomic.fna.gz'

for kingdom in ['archaea', 'bacteria']:
    fna = assembly[assembly['group'] == kingdom]['fna'].to_list()
    n = 8 if kingdom == 'archaea' else 256 # split into 8 or 256 chunks so that each contains same number of genomes
    chunks = [fna[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        with open('assembly/fna/' + kingdom + '.split.' + str(i) + '.id', 'w') as w:
            w.write('\n'.join(chunk) + '\n')

deprecated = (assembly.ftp_path == 'na') | (assembly.ftp_path.isnull())
print(f'#species: {assembly.taxonomy.nunique()}')
print(f'#assemblies: {assembly.assembly.nunique()}')
print(f'#deprecated_reps: {len(assembly[deprecated])}')
"

# ## grab representative genomes if necessary
# wget -qN --show-progress https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
# tar -xvf gtdb_genomes_reps.tar.gz

# python -c "
# import pandas as pd
# import glob
# import os
# import shutil
# dset = set(pd.read_table('assembly/deprecated.tsv').assembly)
# files = glob.glob('gtdb_genomes_reps/**/*.fna.gz', recursive=True)
# files = [file for file in files if os.path.basename(file).rsplit('_genomic')[0] in dset]
# os.makedirs('assembly/fna/deprecated_reps', exist_ok=True)
# for file in files:
#     shutil.copy(file, f'assembly/fna/deprecated_reps/{os.path.basename(file)}')
# "
# 
# cat assembly/fna/deprecated_reps/*.fna.gz > assembly/fna/deprecated_reps.fna.gz
```

### Step 3: Download assemblies and generate an accession2assembly mapping

```bash
## download assemblies then cat
find assembly/fna -maxdepth 1 -name '*.id' | xargs -P 32 -I {} bash -c '
    wget -i ${1} -qN --show-progress -P ${1%.id}; \
    find ${1%.id} -maxdepth 1 -name "*.fna.gz" | xargs cat > ${1%.id}.fna.gz' - {}

## create a dictionary that maps sequence to assembly
python -c "
import glob
import os
import gzip
import pandas as pd
from tqdm.contrib.concurrent import process_map

def parse_file(file):
    lines = []
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line[0] == '>':
                lines.append([line[1:].split()[0], '_'.join(os.path.basename(file).split('_')[:2])])
    return lines

pd.DataFrame([
    x for y in process_map(parse_file, glob.glob('assembly/fna/**/*.fna.gz'), max_workers=64, chunksize=1) for x in y
], columns=['accession', 'assembly']).to_csv('assembly/accession2assembly.tsv', sep='\t', index=False, header=None)
"
```

### Step 4: Collect NCBI plasmids


### Step 5: Collect SARG+
```bash
mkdir -p prot
git clone https://github.com/xinehc/sarg-curation
mv sarg-curation/sarg.fa prot/prot.fa
rm -rf sarg-curation
```

### Step 5 (alternative): Collect NDARO
Protein sequences and metadata of NDARO need to be manually downloaded from https://www.ncbi.nlm.nih.gov/pathogens/refgene/ (click `Download` for both the metadata `refgenes.tsv` and reference protein sequences `refgene_catalog.zip`). After unzip, place all these two files to a folder called `prot`.

We need to parse the metadata to make sure it have a two-level hierarchy.
```bash
```

### Step 5 (alternative): Collect CARD
CARD can be obtained from https://card.mcmaster.ca/download. We use the `protein_fasta_protein_homolog_model.fasta` and the `aro_index.tsv` file only. Please place them to a folder called `prot`.



## Construction of the taxonomic database
### Step 1: Map assemblies to the protein databases
> [!NOTE]
> This step will take a while to finish. If you are working on HPC, please submit a SLURM job for each `*.fna.gz` file and increase `--threads` to reduce computational time. If not, combining all `*.fna.gz` files into a single one then running `diamond` with tuned `-b -c` may help to speed up.

```bash
mkdir -p nucl/out

find assembly/fna -maxdepth 1 -name '*.fna.gz' | sort | xargs -P 8 -I {} bash -c '
    filename=${1%.fna*};
    filename=${filename##*/};
    echo $filename;
    diamond blastx \
        --db prot/prot.fa \
        --query $1 \
        --out nucl/out/${filename}.txt \
        --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
        --evalue 1e-15 --subject-cover 90 --id 90 \
        --range-culling --frameshift 15 --range-cover 25 \
        --max-hsps 0 --max-target-seqs 25 \
        --threads 8 --quiet' - {}
```

### Step 2: Parse output files and extract sequences
Parse output files of `diamond` by ensuring at most 25% overlap of HSPs, keep only a single sequence per qseqid-gene combination, discard also sequences close to the boundaries.

```bash
mkdir -p nucl/seq

python -c "
import os
import glob
import subprocess
import pandas as pd
from math import floor, ceil
from collections import defaultdict
from tqdm.contrib.concurrent import process_map

def sort_coordinate(start, end):
    return (start - 1, end, '+') if start < end else (end - 1, start, '-')

def compute_overlap(coordinates):
    qstart, qend, sstart, send = coordinates
    overlap = min(qend, send) - max(qstart, sstart)
    return max(overlap / (qend - qstart), overlap / (send - sstart))

def extract_sequence(file):
    filename = os.path.basename(file).split('.fna.gz')[0]
    qrange = defaultdict(set)
    lines = []

    with open('nucl/out/' + filename + '.txt') as f:
        for line in f:
            ls = line.split()
            qseqid, sseqid = ls[0], ls[1]
            qstart, qend, strand = sort_coordinate(int(ls[5]), int(ls[6]))
            if (
                qseqid not in qrange or
                all([compute_overlap((qstart, qend, *x)) < 0.25 for x in qrange.get(qseqid)])
            ):
                qrange[qseqid].add((qstart, qend))
                pident = float(ls[2])
                qlen, slen = int(ls[4]), int(ls[7]) * 3
                qcoord = floor((qstart + qend) / 2) if strand == '+' else ceil((qstart + qend) / 2)

                ## make sure not too close to the boundary
                if qend + 2500 < qlen and qstart - 2500 > 0:
                    lines.append([qseqid, qcoord - 5000, qcoord + 5000, sseqid, pident, strand, sseqid.split('-')[0]])

    ## keep only one sequence per qseqid + gene
    lines = pd.DataFrame(lines, columns=['qseqid', 'qstart', 'qend', 'sseqid', 'pident', 'strand', 'gene'])
    lines = lines.sort_values(['qseqid', 'gene', 'pident', 'strand'], ascending=False).groupby(['qseqid', 'gene'], as_index=False).first()
    lines.drop('gene', axis=1).to_csv('nucl/seq/' + filename + '.bed', sep='\t', header=None, index=False)

    ## extract sequneces from original files
    with open('nucl/seq/' + filename + '.fa', 'w') as f:
        subprocess.run([
            'seqkit', 'subseq', file,
            '--bed', 'nucl/seq/' + filename + '.bed'
        ], stdout=f, stderr=subprocess.DEVNULL, check=True)

process_map(extract_sequence, glob.glob('assembly/fna/*.fna.gz'), max_workers=64, chunksize=1)
"
```

### Step 3: Cluster to remove duplicated sequences
Create a file for each species-gene combo.

```bash
mkdir -p nucl/raw
find nucl/seq -maxdepth 1 -name '*.fa' | sort | xargs cat > nucl/raw.fa

python -c "
import pandas as pd
from collections import defaultdict

archaea = {'l2', 'l11', 'l10e', 'l15e', 'l18e', 's3ae', 's19e', 's28e'}
bacteria = {'l2', 'l11', 'l20', 'l27', 's2', 's7', 's9', 's16'}

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

sequence = defaultdict(list)
accession = set()
with open('nucl/raw.fa') as f:
    for line in f:
        if line[0] == '>':
            ls = line.split()
            gs = ls[1].split('-')
            ac = ls[0][1:].rsplit('_', 1)[0]
            species = assembly2species.get(accession2assembly.get(ac))
            phylum = species[0].split(';')[0].split('|')[-1]
            if (
                phylum == 'Archaea' and gs[4] == 'archaea' and gs[0] in archaea or
                phylum == 'Bacteria' and gs[4] == 'bacteria' and gs[0] in bacteria
            ):
                save = True
                accession.add(ac)
                subset = gs[0].replace('/', '_') + '.' + species[-1]
            else:
                save = False
        if save:
            sequence[subset].append(line)

## split the sequences into two parts
with open('nucl/nucl_a.fa', 'w') as w:
    for key, val in sequence.items():
        record = ''.join(val)
        if record.count('>') == 1:
            w.write(record)
        else:
            with open('nucl/raw/{}.fa'.format(key), 'w') as ww:
                ww.write(record)

## save the metadata file for later usage
pd.DataFrame([
    [x] + assembly2species.get(accession2assembly.get(x))[0].split(';') for x in accession
], columns=['accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']).sort_values('accession').to_csv('nucl/raw.tsv', index=False, sep='\t')
"
```

> [!NOTE]
> If your are working on HPC, please consider submitting a SLURM job for each `mmseqs` to make them run in parallel, and increase `--threads` to reduce the computational time. Lower `-P` if you encounter memory issues.

```bash
mkdir -p nucl/clustered

find nucl/raw -maxdepth 1 -name '*.fa' | sort | xargs -P 16 -I {} bash -c '
    filename=${1%.fa*};
    filename=${filename##*/};
    FILE_SIZE=$(stat -c "%s" nucl/raw/$filename.fa)
    if [ "$FILE_SIZE" -gt 10000000 ]; then THREADS=4; else THREADS=1; fi;
    echo $filename $FILE_SIZE $THREADS;
    mmseqs easy-cluster \
        $1 nucl/clustered/$filename nucl/clustered/$filename \
        -c 0.9995 --min-seq-id 0.9995 --cov-mode 1 \
        -s 7.5 --cluster-reassign --threads $THREADS -v 0 > /dev/null' - {}
```

Combine all clustered files to get the nucleotide database.

```bash
find nucl/clustered -maxdepth 1 -name '*rep_seq.fasta' -exec cat {} \; > nucl/nucl_b.fa
cat nucl/nucl_a.fa nucl/nucl_b.fa | seqkit sort | seqkit shuffle -s 0 > nucl/nucl.fa

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
            subset = ls[1].split('-')[-2] + '.' + ls[1].split('-')[0].replace('/', '_')
        sequence[subset].append(line)

for key, val in sequence.items():
    with open('nucl/nucl.{}.fa'.format(key), 'w') as w:
        w.write(''.join(val))

metadata = pd.read_table('nucl/raw.tsv')
metadata[metadata.accession.isin(accession)].to_csv('nucl/metadata.tsv', index=False, sep='\t')
"
```

## Construction of the plasmid database
### Step 0: selection
### Step 1: geNomad prescreening

### Step 2: geNomad postscreening




## Compress
Get all necessary files into the database.

```bash
mkdir -p database
cp nucl/metadata.tsv nucl/nucl.bacteria*.fa nucl/nucl.archaea*.fa database
cat plus/plus_rep_seq.fasta prot/prot.fa | seqkit sort | seqkit shuffle -s 0 > database/prot.fa
tar --sort=name -zcvf database.tar.gz database
```
# argo-supplementary
