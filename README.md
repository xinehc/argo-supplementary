# Argo-supplementary

Instructions for building/extending the database of [Argo](https://github.com/xinehc/argo).

- [Argo-supplementary](#argo-supplementary)
  - [Prerequisites](#prerequisites)
    - [Step 1: Install necessary packages](#step-1-install-necessary-packages)
    - [Step 2: Download GTDB assemblies](#step-2-download-gtdb-assemblies)
    - [Step 3: Download RefSeq *complete genome* and *chromosome-level* plasmids](#step-3-download-refseq-complete-genome-and-chromosome-level-plasmids)
    - [Step 4: Collect SARG+](#step-4-collect-sarg)
    - [(alternative) Step 4: Collect NDARO](#alternative-step-4-collect-ndaro)
    - [(alternative) Step 4: Collect CARD](#alternative-step-4-collect-card)
  - [Construction of the *reference taxonomy* database](#construction-of-the-reference-taxonomy-database)
    - [Step 1: Create an accession to assembly mapping](#step-1-create-an-accession-to-assembly-mapping)
    - [Step 2: Annotate ARGs with `DIAMOND`](#step-2-annotate-args-with-diamond)
    - [Step 3: Parse annotations and extract sequences](#step-3-parse-annotations-and-extract-sequences)
  - [Construction of the *reference plasmid* database](#construction-of-the-reference-plasmid-database)
    - [Step 1: Retrieve putative plasmid sequences from assemblies](#step-1-retrieve-putative-plasmid-sequences-from-assemblies)
    - [Step 2: Extract ARG-containing sequences](#step-2-extract-arg-containing-sequences)
    - [Step 3: Filter out non-plasmid sequences with `geNomad`](#step-3-filter-out-non-plasmid-sequences-with-genomad)
  - [Deduplication](#deduplication)
    - [Step 1: Group ARG-containing sequences for each species](#step-1-group-arg-containing-sequences-for-each-species)
    - [Step 2: Cluster of sequences with `MMseqs2`](#step-2-cluster-of-sequences-with-mmseqs2)
    - [Step 3: Split files according to ARG type](#step-3-split-files-according-to-arg-type)
  - [Compress](#compress)
  - [(optional) Database extension](#optional-database-extension)
    - [Step 1: Install necessary packages and prepare reference genomes](#step-1-install-necessary-packages-and-prepare-reference-genomes)
    - [Step 2: Annotate ARGs](#step-2-annotate-args)
    - [Step 3: Predict taxonomy and genomic context](#step-3-predict-taxonomy-and-genomic-context)
    - [Step 4: Merge sequences and metadata](#step-4-merge-sequences-and-metadata)
## Prerequisites

### Step 1: Install necessary packages

```bash
conda install -c bioconda -c conda-forge 'seqkit>=2.6.1' 'genomad>=1.8.1' 'mmseqs2>=15.6f452' 'diamond==2.1.8' 'tqdm' 'pandas' 'biopython'
genomad download-database .
```

### Step 2: Download GTDB assemblies

> [!NOTE]
> Some genomes may be deprecated by NCBI at the time of GTDB releases. Representative genomes can still be downloaded from the GTDB FTP (`genomic_files_reps/gtdb_genomes_reps.tar.gz`) and included manually, but the non-representative ones can no longer be retrieved.

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

## read metadata, split assemblies into chunks according to phyla
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
```

Download assemblies for each `*.id` from NCBI FTP, then concatenate them into chunks.

```bash
find assembly/fna -maxdepth 1 -name '*.id' | xargs -P 32 -I {} bash -c '
    wget -i ${1} -qN --show-progress -P ${1%.id}; \
    find ${1%.id} -maxdepth 1 -name "*.fna.gz" | xargs cat > ${1%.id}.fna.gz' - {}
```

Grab representative genomes if necessary.

```bash
wget -qN --show-progress https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz
tar -xvf gtdb_genomes_reps.tar.gz

python -c "
import pandas as pd
import glob
import os
import shutil

dset = set(pd.read_table('assembly/deprecated.tsv').assembly)
files = glob.glob('gtdb_genomes_reps/**/*.fna.gz', recursive=True)
files = [file for file in files if os.path.basename(file).rsplit('_genomic')[0] in dset]
os.makedirs('assembly/fna/deprecated_reps', exist_ok=True)
for file in files:
    shutil.copy(file, f'assembly/fna/deprecated_reps/{os.path.basename(file)}')
"

cat assembly/fna/deprecated_reps/*.fna.gz > assembly/fna/deprecated_reps.fna.gz
```

### Step 3: Download RefSeq *complete genome* and *chromosome-level* plasmids

```bash
mkdir -p plasmid/fna
wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -P plasmid

python -c "
import pandas as pd

assembly = pd.read_table('assembly/assembly_summary_refseq.txt', skiprows=1, low_memory=False).rename({'#assembly_accession': 'assembly'}, axis=1)
assembly = assembly[(assembly.group.isin(['archaea', 'bacteria'])) & ((assembly['ftp_path'] != 'na'))]
assembly = assembly[assembly.assembly_level.isin({'Complete Genome', 'Chromosome'})]
assembly['fna'] = assembly['ftp_path'] + '/' + assembly['ftp_path'].str.split('/').str.get(-1) + '_genomic.fna.gz'

for kingdom in ['archaea', 'bacteria']:
    fna = assembly[assembly['group'] == kingdom]['fna'].to_list()
    n = 2 if kingdom == 'archaea' else 32
    chunks = [fna[i::n] for i in range(n)]

    for i, chunk in enumerate(chunks):
        with open('plasmid/fna/' + kingdom + '.split.' + str(i) + '.id', 'w') as w:
            w.write('\n'.join(chunk) + '\n')

print(f'#assemblies: {assembly.assembly.nunique()}')
"
```

Download assemblies.

```bash
find plasmid/fna -maxdepth 1 -name '*.id' | xargs -P 32 -I {} bash -c '
    wget -i ${1} -qN --show-progress -P ${1%.id}; \
    find ${1%.id} -maxdepth 1 -name "*.fna.gz" | xargs cat > ${1%.id}.fna.gz' - {}
```

### Step 4: Collect SARG+

> [!NOTE]
> Argo uses [SARG+](https://github.com/xinehc/sarg-curation) by default for ARG annotation, but it also supports [NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/) (Collect National Database of Antibiotic Resistant Organisms) and [CARD](https://card.mcmaster.ca/) (Comprehensive Antibiotic Resistance Database) (see details below).

```bash
mkdir -p prot
git clone https://github.com/xinehc/sarg-curation
mv sarg-curation/sarg.fa prot/prot.fa
```

### (alternative) Step 4: Collect NDARO

Reference protein sequences and metadata of NDARO can be downloaded from https://www.ncbi.nlm.nih.gov/pathogens/refgene/ (click `Download` for both the metadata `refgenes.tsv` and the reference protein sequences `protein.faa`). After unzipping, place both files into a folder named `prot`.

We need to parse the metadata to ensure it has a two-level hierarchy.

```bash
python -c "
import pandas as pd
from Bio import SeqIO

ndaro = pd.read_table('prot/refgenes.tsv').fillna('NA')
ndaro['accession'] = ndaro.apply(lambda x: x['RefSeq protein'] if x['RefSeq protein']!='NA' else x['GenBank protein'], axis=1)
ndaro = ndaro[(ndaro.groupby('accession').transform('size') == 1) & (ndaro.Subtype.isin(['BIOCIDE', 'AMR']))]
ndaro['id'] = 'NDARO|' + ndaro['Class'].str.lower().str.replace(' ', '_') + '|' + ndaro['Gene family'].str.replace(' ', '_')
acc2id = ndaro.set_index('accession')['id'].to_dict()

records = []
dups = set()
with open('prot/proteins.faa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        if record.id in acc2id:
            record.id = acc2id.get(record.id) + '|' + record.id
            if record.seq not in dups:
                dups.add(record.seq)
                records.append(record)

with open('prot/prot.fa', 'w') as output_handle:
    SeqIO.write(records, output_handle, 'fasta')
"
```

### (alternative) Step 4: Collect CARD

CARD can be obtained from https://card.mcmaster.ca/download. Files `protein_fasta_protein_homolog_model.fasta` and `aro_index.tsv` should be placed in a folder named `prot`.

We need to parse the metadata to ensure it has a two-level hierarchy.

```bash
python -c "
import pandas as pd
import re
from Bio import SeqIO

card = pd.read_table('prot/aro_index.tsv').fillna('NA')
card['id'] = 'CARD|' + (card['Drug Class'].str.lower().str.replace(' antibiotic', '') + '|' + card['CARD Short Name']).str.replace(' ', '_')
acc2id = card.set_index('ARO Accession')['id'].to_dict()

records = []
dups = set()
with open('prot/protein_fasta_protein_homolog_model.fasta') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        record.id = acc2id.get(re.search('ARO:[0-9]+', record.description).group()) + '|' + record.id.split('|')[-3]
        if record.id and record.seq not in dups:
            dups.add(record.seq)
            records.append(record)

with open('prot/prot.fa', 'w') as output_handle:
    SeqIO.write(records, output_handle, 'fasta')
"
```

## Construction of the *reference taxonomy* database

### Step 1: Create an accession to assembly mapping

Create a dictionary that maps sequence accessions to to their corresponding assemblies.

```bash
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

### Step 2: Annotate ARGs with `DIAMOND`

> [!NOTE]
> This step will take a while to finish. If you are working on an HPC, please submit a separate SLURM job for each `*.fna.gz` file and increase `--threads` to reduce computation time. If not, combining all `*.fna.gz` files into a single one and running `diamond` with tuned `-b -c` may help speed up the process.

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

### Step 3: Parse annotations and extract sequences

Parse output files of `diamond`, ensuring at most 25% pairwise overlap of ARGs. Discard also ARGs close to the boundaries, except for sequences with circular topology.

```bash
mkdir -p nucl/seq 

python -c "
import pandas as pd
import glob
import gzip
import os
import xml.etree.ElementTree as ET

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

## retrieve sequence topology from NCBI
def get_topology(ids, chunksize=1000):
    topology = []
    n = len(ids) // chunksize + 1
    folds = [list(ids)[i::n] for i in range(n)]
    for fold in tqdm(folds):
        tmp = []
        while len(tmp) != len(fold):
            with Entrez.esummary(db='nuccore', id=','.join(fold), rettype='xml', version='2.0') as handle:
                root = ET.parse(handle).getroot()
                tmp = [[i,j] for i,j in zip([x.text for x in root.findall('.//AccessionVersion')], [x.text for x in root.findall('.//Topology')])]

            try: assert len(tmp) == len(fold)
            except: 
                ## much slower with efetch but safer in general
                todo = set(fold) - set(x[0] for x in tmp)
                with Entrez.efetch(db='nuccore', id=','.join(todo), rettype='gb', retmode='xml') as handle:
                    for i in Entrez.read(handle):
                        tmp.append([i['GBSeq_accession-version'], i['GBSeq_topology']])
        topology.extend(tmp)

    return topology

## extract sequence from GTDB
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

                ## make sure not too close to the boundary
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

## get ids of complete genome and chromosome-level assemblies
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

## record sequences with ARGs close to the boundary
ids = set()
for file in tqdm(glob.glob('nucl/out/*.txt')):
    with open(file) as f:
        for line in f:
            ls = line.rstrip().split()
            qlen, qstart, qend = int(ls[4]), int(ls[5]), int(ls[6])
            if max(qstart, qend) + 5000 > qlen or min(qstart, qend) - 5000 < 0:
                if ls[0] in aset:
                    ids.add(ls[0])

qseqid2topology = {x[0]: x[1] for x in get_topology(ids)}
process_map(extract_sequence, glob.glob('assembly/fna/*.fna.gz'), max_workers=64, chunksize=1)
"
```

## Construction of the *reference plasmid* database

### Step 1: Retrieve putative plasmid sequences from assemblies

Extract sequences containing keywords `plasmid` or `megaplasmid`.

```bash
python -c "
import glob
import gzip
from Bio import SeqIO
from tqdm.contrib.concurrent import process_map

def parser(file):
    records = []
    with gzip.open(file, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if ' plasmid' in record.description or ' megaplasmid' in record.description:
                records.append(record)
    return records

r = process_map(parser, glob.glob('plasmid/fna/*/*.fna.gz'), max_workers=64, chunksize=1)
with open('plasmid/plasmid.fna', 'w') as output_handle:
    SeqIO.write([x for y in r for x in y], output_handle, 'fasta')
"
```

### Step 2: Extract ARG-containing sequences

```bash
diamond blastx \
    --db prot/prot.fa \
    --query plasmid/plasmid.fna \
    --out plasmid/plasmid_sarg.txt \
    --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
    --evalue 1e-15 --subject-cover 90 --id 90 \
    --range-culling --frameshift 15 --range-cover 25 \
    --max-hsps 0 --max-target-seqs 25 \
    --threads 64 --quiet
```

Extract ARG-containing sequences based on annotated coordinates.

```bash
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
    for fold in tqdm(folds):
        tmp = []
        while len(tmp) != len(fold):
            with Entrez.esummary(db='nuccore', id=','.join(fold), rettype='xml', version='2.0') as handle:
                root = ET.parse(handle).getroot()
                tmp = [[i,j] for i,j in zip([x.text for x in root.findall('.//AccessionVersion')], [x.text for x in root.findall('.//Topology')])]

            try: assert len(tmp) == len(fold)
            except: 
                ## much slower with efetch but safer in general
                todo = set(fold) - set(x[0] for x in tmp)
                with Entrez.efetch(db='nuccore', id=','.join(todo), rettype='gb', retmode='xml') as handle:
                    for i in Entrez.read(handle):
                        tmp.append([i['GBSeq_accession-version'], i['GBSeq_topology']])
        topology.extend(tmp)

    return topology

## record sequences with ARGs close to the boundary
ids = set()
with open('plasmid/plasmid_sarg.txt') as f:
    for line in f:
        ls = line.rstrip().split()
        qlen, qstart, qend = int(ls[4]), int(ls[5]), int(ls[6])
        if max(qstart, qend) + 5000 > qlen or min(qstart, qend) - 5000 < 0:
            ids.add(ls[0])

## retrieve sequence topology from NCBI
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

            ## make sure not too close to the boundary
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
```


### Step 3: Filter out non-plasmid sequences with `geNomad`

RefSeq plasmids can be mislabeled, use `geNomad` for double-checking.

```bash
## first round for filtering chromosomes
seqkit grep -f <(cut plasmid/plasmid_sarg.txt -f1) plasmid/plasmid.fna > plasmid/plasmid_sub.fa
genomad end-to-end --cleanup plasmid/plasmid_sub.fa plasmid/genomad genomad_db \
    --disable-find-proviruses \
    --conservative \
    --threads 64

## second round with permisssive cutoffs for filtering chimeras
genomad end-to-end --cleanup plasmid/plasmid_arg.fa plasmid/genomad genomad_db \
    --disable-find-proviruses \
    --relaxed \
    --threads 64
```

Obtain a subset of ARG-containing plasmid sequences based on `geNomad`'s results.

```bash
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
```

## Deduplication

### Step 1: Group ARG-containing sequences for each species

Merge all files into a single one.

```bash
mkdir -p nucl/raw
find nucl/seq -maxdepth 1 -name '*.fa' | sort | xargs cat > nucl/raw.fa
```

Create a file for each species-type-subtype combo.

```bash
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
            subset = ls[1].split('|')[1] + '@' + ls[1].split('|')[2].replace('/','_SLASH_').replace('*','_STAR').replace('''\'''', '_PRIME')+ '.' + lineage
        sequence[subset].append(line)

## split the sequences into two parts
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

## save the metadata file for later usage
with open('nucl/raw.tsv', 'w') as f:
    f.write('\t'.join(['accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) + '\n')
    for i in accession:
        f.write('\t'.join([i] + assembly2species.get(accession2assembly.get(i))[0].split(';')) + '\n')
"
```

### Step 2: Cluster of sequences with `MMseqs2`

> [!NOTE]
> This step will take some time to complete. If you are working on an HPC, you may consider splitting `nucl/raw/*.fa` into chunks and submit a SLURM job for each chunk to run them in parallel. Increase `--threads` to reduce computational time, and lower `-P` if you encounter memory issues.

```bash
mkdir -p nucl/clustered

find nucl/raw -maxdepth 1 -name '*.fa' | sort | xargs -P 64 -I {} bash -c '
    filename=${1%.fa*}
    filename=${filename##*/}
    FILE_SIZE=$(stat -c "%s" $1)

    if [ ! -f nucl/clustered/${filename}_rep_seq.fasta ]; then
        if [ $FILE_SIZE -gt 10000000 ]; then THREADS=16; else THREADS=1; fi
        echo $filename $FILE_SIZE $THREADS;
        mmseqs easy-cluster \
            $1 nucl/clustered/$filename nucl/clustered/$filename \
            -c 0.9995 --min-seq-id 0.9995 --cov-mode 1 \
            -s 7.5 --cluster-reassign --threads $THREADS -v 0 > /dev/null
        rm -rf nucl/clustered/$filename nucl/clustered/${filename}_all_seqs.fasta nucl/clustered/${filename}_cluster.tsv
    fi' - {}

python -c "
import glob

a = {x.split('/')[-1].split('_rep_seq.fasta')[0] for x in glob.glob('nucl/clustered/*_rep_seq.fasta')}
b = {x.split('/')[-1].split('.fa')[0] for x in glob.glob('nucl/raw/*.fa')}

print(len(a), len(b))
assert a==b, 'Some files are not processed for some reasons. Try again.'
"
```

### Step 3: Split files according to ARG type

Combine all clustered files to get the database.

```bash
find nucl/clustered -maxdepth 1 -name '*rep_seq.fasta' -exec cat {} \; > nucl/nucl_c.fa
cat nucl/nucl_*.fa | seqkit sort | seqkit shuffle -s 0 > nucl/nucl.fa

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
        sequence[subset].append(line)

for key, val in sequence.items():
    with open('nucl/sarg.{}.fa'.format(key), 'w') as w:
        w.write(''.join(val))

metadata = pd.read_table('nucl/raw.tsv')
metadata[metadata.accession.isin(accession)].sort_values('accession').to_csv('nucl/sarg.metadata.tsv', index=False, sep='\t')
"
```

## Compress

Get all necessary files into the database.

```bash
wget -qN --show-progress https://zenodo.org/records/12571554/files/database.tar.gz
tar -zxvf database.tar.gz

cp prot/prot.fa database/sarg.fa
cp nucl/sarg.metadata.tsv nucl/sarg*.fa database

tar --sort=name -zcvf database.tar.gz database
```

## (optional) Database extension

### Step 1: Install necessary packages and prepare reference genomes

To incorporate additional sequences into the reference taxonomy/plasmid databases, we need `GTDB-Tk` for taxonomic classification.

```bash
conda install -c bioconda -c conda-forge 'gtdbtk'
download-db.sh gtdbtk_db
```

Using the assemblies of a wildtype *Escherichia coli* and a wildtype *Enterococcus lactis* (isolated from HK WWTP) as an example:

```bash
wget -qN --show-progress https://zenodo.org/records/13992057/files/genome.zip -P extension
unzip extension/genome.zip -d extension

## rename to ensure unique contig ID
python -c "
import glob
import os
from Bio import SeqIO

records = []
topology = []
for file in glob.glob('extension/genome/*.fasta'):
    newfile = file.replace(' ', '_')
    os.rename(file, newfile)
    with open(newfile) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            record.id = file.split('/')[-1].split('.fasta')[0].replace(' ', '_') + '|' + record.id
            records.append(record)

            ## if topology is unknown, simply set topology to 'linear'
            topology.append([record.id, 'linear' if 'circular=true' not in record.description else 'circular'])
        
with open('extension/genome.fna', 'w') as w:
    SeqIO.write(records, w, 'fasta')

with open('extension/topology.txt', 'w') as w:
    for row in topology:
        w.write('\t'.join(row) + '\n')
"
```

### Step 2: Annotate ARGs

```bash
diamond blastx \
    --db prot/prot.fa \
    --query extension/genome.fna \
    --out extension/genome_sarg.txt \
    --outfmt 6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore \
    --evalue 1e-15 --subject-cover 90 --id 90 \
    --range-culling --frameshift 15 --range-cover 25 \
    --max-hsps 0 --max-target-seqs 25 \
    --threads 64 --quiet
```

```bash
python -c "
from collections import defaultdict
from math import floor, ceil
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def sort_coordinate(start, end):
    return (start - 1, end, '+') if start < end else (end - 1, start, '-')

def compute_overlap(coordinates):
    qstart, qend, sstart, send = coordinates
    overlap = min(qend, send) - max(qstart, sstart)
    return max(overlap / (qend - qstart), overlap / (send - sstart))

qseqid2topology = dict()
with open('extension/topology.txt') as f:
    for line in f:
        ls = line.rstrip().split('\t')
        qseqid2topology[ls[0]] = ls[1]

bed = defaultdict(list)
qrange = defaultdict(set)
with open('extension/genome_sarg.txt') as f:
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

            ## make sure not too close to the boundary
            topology = qseqid2topology.get(qseqid)
            if (qend + 2500 < qlen and qstart - 2500 > 0) or topology == 'circular':
                bed[ls[0]].append((qcoord - 5000, qcoord + 5000, strand, sseqid, topology))

records = []
with open('extension/genome.fna') as handle:
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

with open('extension/genome_arg.fa', 'w') as output_handle:
    SeqIO.write(records, output_handle, 'fasta')
"
```

### Step 3: Predict taxonomy and plasmids

```bash
gtdbtk classify_wf --genome_dir extension/genome/ --cpus 64 --out_dir extension/gtdbtk -x fasta --mash_db gtdbtk_db/release220/mash/
seqkit grep -f <(cut extension/genome_sarg.txt -f1) extension/genome.fna > extension/genome_sub.fa
genomad end-to-end --cleanup extension/genome_sub.fa extension/genomad genomad_db \
    --disable-find-proviruses \
    --conservative \
    --threads 64

## second round with permisssive cutoffs for filtering chimeras
genomad end-to-end --cleanup extension/genome_arg.fa extension/genomad genomad_db \
    --disable-find-proviruses \
    --relaxed \
    --threads 64
```

```bash
python -c "
import pandas as pd
from Bio import SeqIO

fea = pd.read_table('extension/genomad/genome_sub_marker_classification/genome_sub_features.tsv')
sub = pd.read_table('extension/genomad/genome_sub_summary/genome_sub_plasmid_summary.tsv')
sub = set(sub[sub.seq_name.isin(fea[fea.n_uscg == 0].seq_name)].seq_name)

arg = pd.read_table('extension/genomad/genome_arg_aggregated_classification/genome_arg_aggregated_classification.tsv')
arg = arg[arg.seq_name.str.split('@').str.get(-1).str.rsplit('_', n=1).str.get(0).isin(sub)]
arg = arg[arg.plasmid_score / arg.chromosome_score > 2]
arg = set(arg.seq_name)

records = []
with open('extension/genome_arg.fa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        if record.id in arg:
            record.id = 'plasmid@' + record.id
            record.description = record.description.split(' ', 1)[-1]
            records.append(record)

with open('extension/plasmid.fa', 'w') as output_handle:
    SeqIO.write(records, output_handle, 'fasta')
"
```

### Step 4: Merge sequences and metadata

To ensure consistency, all genomes that lack *species-level* resolution will not be used.

```bash
python -c "
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import os

files = []
for file in ['gtdbtk.ar53.summary.tsv', 'gtdbtk.bac120.summary.tsv']:
    file = 'extension/gtdbtk/' + file
    if os.path.isfile(file):
        files.append(pd.read_table(file))

gtdbtk = pd.concat(files)
gtdbtk = gtdbtk[gtdbtk.classification.str.split(';s__').str.get(-1) != '']
gtdbtk['taxonomy'] = gtdbtk.classification.str.replace('[a-z]__', '', regex=True)
filename2taxonomy = gtdbtk.set_index('user_genome').taxonomy.to_dict()

meta = set()
records = defaultdict(list)
with open('extension/genome_arg.fa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        filename = record.id.split('|')[0]
        meta.add((record.id.rsplit('_', 1)[0], *filename2taxonomy.get(filename).split(';')))
        records[record.description.split(' ')[1].split('|')[1]].append(record)

with open('extension/plasmid.fa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        records[record.description.split(' ')[1].split('|')[1]].append(record)

for i, j in records.items():
    with open(f'extension/sarg.{i}.fa', 'w') as w:
        SeqIO.write(j, w, 'fasta')

with open('extension/sarg.metadata.tsv', 'w') as w:
    for line in meta:
        w.write('\t'.join(line) + '\n')
"
```

Merge newly extracted sequences with existing databases.

```bash
mkdir -p database.extension
for file in database/sarg.*.fa database/sarg.metadata.tsv
do
    filename=${file##*/}
    echo $filename
    if [ ! -e extension/$file ]; then
        cp database/$filename database.extension/$filename
    else 
        cat database/$filename extension/$filename > database.extension/$filename
    fi
done
```
