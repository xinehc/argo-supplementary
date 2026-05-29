# Argo-supplementary

Instructions for building the database of [Argo](https://github.com/xinehc/argo).

#### Install necessary packages
```bash
conda install -c bioconda -c conda-forge 'python=3.11' 'seqkit>=2.13.0' 'genomad>=1.12.0' 'mmseqs2>=18.8cc5c' 'diamond>=2.2.0' 'tqdm' 'pandas' 'biopython' 'tar' 'wget'
genomad download-database .
```

#### Build the GTDB database
```bash
bash run.sh
```