#!/usr/bin/env bash

cd data/encode_yy1_zhx

xargs -L 1 curl -O -J -L < encode_manifest.txt

cd ../ENCODE_YY1_293T

xargs -L 1 curl -O -J -L < files.txt