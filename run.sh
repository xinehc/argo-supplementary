#!/usr/bin/env bash
set -euo pipefail

## download necessary files
bash cmd/1_download.sh
bash cmd/2_download.sh
bash cmd/3_download.sh

## build GTDB database
bash cmd/4_annotate.sh
bash cmd/5_annotate.sh

## cluster and compress
bash cmd/6_cluster.sh
bash cmd/7_compress.sh
