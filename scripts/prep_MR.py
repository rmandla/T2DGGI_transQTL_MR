import pandas as pd
import sys

def prep_GWAS_data(gwas_path,protname,dataset,clumped_snps=None,output_suffix=''):
    metal = pd.read_table(gwas_path)

    metal = metal[['CHR','Pos','RSID','Allele1','Allele2','P','Effect','StdErr','Freq1','N','MarkerName_HG38']]
    metal.columns = ['chr','pos','rsid','eAllele','oAllele','P','Effect','SE','AF','N','snpid']
    metal['eAllele'] = metal['eAllele'].str.upper()
    metal['oAllele'] = metal['oAllele'].str.upper()

    if type(clumped_snps) == str:
        snps = pd.read_table(clumped_snps)
        snps['snpid'] = 'chr'+snps['CHR'].astype(str)+':'+snps['BP'].astype(str)
        metal = metal[metal['snpid'].isin(snps['snpid'])]
    metal.drop_duplicates().sort_values('rsid').drop(columns=['snpid']).to_csv(f'{protname}_{dataset}_{output_suffix}_GWAS_harmonized.txt',sep='\t',index=None)

def prep_pQTL_data(protname,dataset,pQTL_path,clumped_snps=None,output_suffix=''):
    gwas = pd.read_table("IGFBP2_decode_GWAS_harmonized.txt")
    pqtl = pd.read_table(pQTL_path)
    
    dataset = dataset.lower()
    if dataset == 'decode':
        pqtl = pqtl[pqtl["rsids"].isin(gwas['rsid'])]
        pqtl['chr'] = pqtl['Chrom'].str.replace('chr','')

        pqtl = pqtl[['chr','Pos','rsids','effectAllele','otherAllele','Pval','Beta','SE','ImpMAF','N']]
    elif dataset == 'ukb':
        gwas['SNP'] = gwas['chr'].astype(str)+':'+gwas['pos'].astype(str)+':'+gwas['eAllele']+':'+gwas['oAllele']
        rsid = gwas[['rsid','SNP']]

        pqtl['SNP1'] = pqtl['CHROM'].astype(str)+':'+pqtl['GENPOS'].astype(str)+':'+pqtl['ALLELE1']+':'+pqtl['ALLELE0']
        pqtl['SNP2'] = pqtl['CHROM'].astype(str)+':'+pqtl['GENPOS'].astype(str)+':'+pqtl['ALLELE0']+':'+pqtl['ALLELE1']

        pqtl1 = pqtl[pqtl['SNP1'].isin(metal['SNP'])]
        pqtl1['SNP'] = pqtl['SNP1']
        pqtl2 = pqtl[pqtl['SNP2'].isin(metal['SNP'])]
        pqtl2['SNP'] = pqtl['SNP2']

        pqtl1 = pqtl1.merge(rsid,left_on='SNP1',right_on='SNP')
        pqtl2 = pqtl2.merge(rsid,left_on='SNP2',right_on='SNP')

        pqtl = pd.concat([pqtl1,pqtl2])
        pqtl['Pval'] = 10**-pqtl['LOG10P']

        pqtl = pqtl[['CHROM','GENPOS','rsid','ALLELE1','ALLELE0','Pval','BETA','SE','A1FREQ','N']]
    else:
        print(f"ERROR: only decode or ukb pQTL dataset names are supported")

    pqtl.columns = ['chr','pos','rsid','eAllele','oAllele','P','Effect','SE','AF','N']
    if type(clumped_snps) == str:
        snps = pd.read_table(clumped_snps)
        snps['snpid'] = snps['CHR'].astype(str)+':'+snps['BP'].astype(str)
        pqtl['snpid'] = pqtl['chr'].astype(str)+':'+pqtl['pos'].astype(str)
        pqtl = pqtl[pqtl['snpid'].isin(snps['snpid'])]
        pqtl = pqtl.drop(columns=['snpid'])
    pqtl.sort_values('rsid').to_csv(f'{protname}_{dataset}_{output_suffix}_pQTL_harmonized.txt',sep='\t',index=None)
