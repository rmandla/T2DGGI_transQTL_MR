# T2DGGI trans Mendelian Randomization Scripts

This repository contains pipelines used to run MR for the T2DGGI GWAS summary statistics against QTL datasets.

## Setup

`git clone https://github.com/rmandla/T2DGGI_transQTL_MR.git`

### Requirements

* Python 3
  * pandas
* R (This code was tested on R 4.1)
  * TwoSampleMR
  * magrittr
  * tidyverse
  * vroom
  * MRPRESSO
* PLINK (This code was tested on PLINK1)

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

**Test for MR in the direction of GWAS -> pQTL**
```
python T2DGGI_transQTL_MR/scripts/prep_MR.py run_mr \
    -s ../MR/EUR_MetalFixed_LDSC-CORR_Results1TBL-HG38-rsid-v2.gz \
    -p IGFBP2 \
    -d decode \
    -ex GWAS \
    -oc pqtl \
    -dpath ../MR/8469_41_IGFBP2_IGFBP_2.txt.gz \
    -cs TEST.GWAS.T2DGGI.ALL.clumped \
    -o MR_output/TEST
```

Here, we use the clumped output above to run MR only using the clumped SNPs. `-s` refers to the input GWAS summary statistic file. `-p` refers to the name of the protein, which is used in naming the output file. `-d` is the pQTL dataset name, which is used to parse the input pQTL dataset. `-ex` sets the exposure type and `-oc` sets the outcome type. These values must be GWAS or PQTL. `-dpath` is the path to the input pQTL summary statistic file. `-cs` is the path to the clumped SNPs generated in the clumping step above. And `-o` is the output header of the MR files. 

This script works by first reformatting and subsetting the input summary statistics to only include the clumped SNPs and all information necessary for MR. Then, it runs an Rscript to get the MR output.
