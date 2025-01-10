# T2DGGI trans Mendelian Randomization Scripts

This repository contains pipelines used to run MR for the T2DGGI GWAS summary statistics against QTL datasets.

## Tutorial

### Genome-wide analyses

#### Data preparation

**Run clumping of imput GWAS summary statistics**
'''
python T2DGGI_transQTL_MR/scripts/prep_MR.py clump \
  -s ../MR/EUR_MetalFixed_LDSC-CORR_Results1TBL-HG38-rsid-v2.gz \
  -r /humgen/florezlab/users/rmandla/1000G \
  -o TEST \
  -e GWAS \
  -p plink
'''

Here, clumping is run against input GWAS in HG38 using reference data from 1000G, where `-r` represents the path containing the reference files in PLINK bim/bed/fam format and `-s` represents the input summary statistics. `-o` represents the output header, such that the output files from the preparation will be `${output_header}.ALL.clumped`. `-e` represents the exposure type, and `-p` is the path to an executable PLINK binary. Since this step runs clumping on the input GWAS, and we are testing for MR between the GWAS and different proteins, this step only needs to be run once for all genome-wide GWAS -> pQTL MR tests.

**Run clumping of input pQTL summary statistics**
``
python T2DGGI_transQTL_MR/scripts/prep_MR.py clump \
  -s ../MR/8469_41_IGFBP2_IGFBP_2.txt.gz \
  -r /humgen/florezlab/users/rmandla/1000G \
  -o TEST2_IGFBP2 \
  -e pQTL \
  -d decode \
  -p plink
``

Here, clumping is run against an input pQTL dataset. The arguments are the same as above, with the `-s` now pointing to the path of the pQTL summary statistics. A new argument, `-d` is added to specify whether the pQTL dataset is deCODE or UKB, which changes how the input `-s` file is parsed.

#### Running MR


