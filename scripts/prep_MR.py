import pandas as pd
import sys

def prep_pQTL_data(protname,dataset,pQTL_path):
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
    pqtl.sort_values('rsid').to_csv(f'{protname}_{dataset}_pQTL_harmonized.txt',sep='\t',index=None)
