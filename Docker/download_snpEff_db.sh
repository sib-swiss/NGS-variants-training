#!/usr/bin/env bash

cd /data
sudo wget https://ngs-variants-training.s3.eu-central-1.amazonaws.com/snpEff_v5_0_GRCh38.99.zip
sudo unzip snpEff_v5_0_GRCh38.99.zip
cd data
sudo mv GRCh38.99/ ../
cd ..
sudo rm snpEff_v5_0_GRCh38.99.zip
sudo rm -r data
