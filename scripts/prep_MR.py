import pandas as pd
import sys, subprocess, argparse, os

def run_clumping(sst,ref_path,exposure,output_header='',dataset=None,plink='plink',snps=None):
    sst_df = pd.read_table(sst)
    if exposure.lower() == 'pqtl':
        if dataset.lower() == 'ukb':
            chrom_col = 'CHROM'
            bp_col = 'GENPOS'
            logp_col = 'LOG10P'
            p_col = None
            a1_col = 'ALLELE0'
            a2_col = 'ALLELE1'
        elif dataset.lower() == 'decode':
            sst_df['CHR'] = sst_df['Chrom'].str.split('chr',expand=True)[1]
            chrom_col = 'CHR'
            bp_col = 'Pos'
            p_col = 'Pval'
            logp_col = None
            a1_col = 'effectAllele'
            a2_col = 'otherAllele'
        else:
            raise(f'ERROR: {dataset} not recognized')
        if output_header != '':
            output_header += '.'
        output = f'{output_header}{exposure}.{dataset}'
    elif exposure.lower() == 'gwas':
        chrom_col = 'CHR'
        bp_col = 'Pos'
        p_col = 'P'
        logp_col = None
        a1_col = 'Allele1'
        a2_col = 'Allele2'
        if output_header != '':
            output_header += '.'
        output = f'{output_header}{exposure}.T2DGGI'
    else:
        raise(f'ERROR: {exposure} not recognized')

    if type(snps) == str:
        snps_df = pd.read_table(snps,header=None)
        sst_df = sst_df[sst_df['MarkerName'].isin(snps_df[0])]

    sst_df['CHR'] = sst_df[chrom_col]
    sst_df = sst_df[sst_df['CHR'].isin([str(i) for i in range(1,23)]+[i for i in range(1,23)])]
    sst_df['CHR'] = sst_df['CHR'].astype(int)
    sst_df['BP'] = sst_df[bp_col]
    if type(logp_col) == str:
        sst_df['P'] = 10**-sst_df[logp_col]
    else:
        sst_df['P'] = sst_df[p_col]
    sst_df['SNP1'] = sst_df['CHR'].astype(str)+':'+sst_df['BP'].astype(str)+':'+sst_df[a1_col].str.upper()+':'+sst_df[a2_col].str.upper()
    sst_df['SNP2'] = sst_df['CHR'].astype(str)+':'+sst_df['BP'].astype(str)+':'+sst_df[a2_col].str.upper()+':'+sst_df[a1_col].str.upper()
    sst_df1 = sst_df.copy()
    sst_df1['SNP'] = sst_df1['SNP1']
    sst_df2 = sst_df.copy()
    sst_df2['SNP'] = sst_df2['SNP2']
    sst_df = pd.concat([sst_df1,sst_df2])
    if '/' in output:
        output_dir = output.sep('/')[:-1].join('/')
        output_file = output.sep('/')[-1]
    else:
        output_dir = './'
        output_file = output
    if output+".ALL.clumped" in os.listdir(output_dir):
        subprocess.run(f'rm {output}.ALL.clumped',shell=True,check=True)
    for chrom in range(1,23):
        tsst_df = sst_df[sst_df['CHR']==chrom]
        tsst_df[['CHR','SNP','BP','P']].to_csv('temp_forclump.txt',sep='\t',index=None)
        ref_file = f'{ref_path}/EUR_1KG_chr{chrom}'
        subprocess.run(f'{plink} --bfile {ref_file} --clump temp_forclump.txt --clump-p1 5e-8 --clump-kb 10000 --clump-r2 0.001 --out {output}.chr{chrom}',shell=True,check=True)

        if f'{output_file}.chr{chrom}.clumped' in os.listdir(output_dir):
            subprocess.run("tail -n+2 "+output+".chr"+str(chrom)+".clumped | awk '{print $3}' | grep -v '^$' >> "+output+".ALL.clumped",shell=True,check=True)

def prep_GWAS_data(gwas_path,protname,dataset,output_header,clumped_snps=None):
    metal = pd.read_table(gwas_path)

    metal = metal[['CHR','Pos','RSID','Allele1','Allele2','P','Effect','StdErr','Freq1','N','MarkerName_HG38']]
    metal.columns = ['chr','pos','rsid','eAllele','oAllele','P','Effect','SE','AF','N','snpid']
    metal['eAllele'] = metal['eAllele'].str.upper()
    metal['oAllele'] = metal['oAllele'].str.upper()

    if type(clumped_snps) == str:
        snps = pd.read_table(clumped_snps)
        snps['snpid'] = 'chr'+snps['CHR'].astype(str)+':'+snps['BP'].astype(str)
        metal = metal[metal['snpid'].isin(snps['snpid'])]
    metal.drop_duplicates().sort_values('rsid').drop(columns=['snpid']).to_csv(f'{output_header}_{protname}_{dataset}_GWAS_harmonized.txt',sep='\t',index=None)

def prep_pQTL_data(protname,dataset,pQTL_path,output_header,clumped_snps=None):
    gwas = pd.read_table(f'{output_header}_{protname}_{dataset}_GWAS_harmonized.txt')
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
    pqtl.sort_values('rsid').to_csv(f'{output_header}_{protname}_{dataset}_pQTL_harmonized.txt',sep='\t',index=None)

def main():
    parser = argparse.ArgumentParser()
    #subparser = parser.add_subparsers()

    subparsers = parser.add_subparsers(title='Commands', dest='command')
    # use dispatch pattern to invoke method with same name

    #parser.add_argument('pull_ld')
    clump = subparsers.add_parser('clump')
    clump.add_argument("-s", "--sst",dest='sst')
    clump.add_argument('-r','--ref_path',dest='ref_path')
    clump.add_argument('-o','--output_header',dest='output_header',default='')
    clump.add_argument('-e','--exposure',dest='exposure')
    clump.add_argument('-d','--dataset',dest='dataset',default=None)
    clump.add_argument('-p','--plink_path',dest='plink_path',default='plink')
    clump.add_argument('-snps','--snps',dest='snps',default=None)

    #parser.add_argument('-c','--compare',dest='compare')
    mr = subparsers.add_parser('run_mr')
    mr.add_argument('-s','--sst',dest='gwas_path')
    mr.add_argument("-p", "--protein_name",dest='protname')
    mr.add_argument('-d', '--dataset',dest='dataset')
    mr.add_argument('-ex','--exposure',dest='exposure')
    mr.add_argument('-oc','--outcome',dest='outcome')
    mr.add_argument('-dpath','--dataset_path',dest='pQTL_path')
    mr.add_argument('-cs','--clumped_snps',dest='clumped_snps',default=None)
    mr.add_argument("-o", "--output",dest = "output_header")

    args = parser.parse_args()
    if args.command == 'clump':
        print('Clumping variants')
        run_clumping(sst=args.sst,ref_path=args.ref_path,exposure=args.exposure,output_header=args.output_header,dataset=args.dataset,plink=args.plink_path,snps=args.snps)
    elif args.command == 'run_mr':
        print('Running MR')
        prep_GWAS_data(gwas_path=args.gwas_path,protname=args.protname,dataset=args.dataset,output_header=args.output_header,clumped_snps=args.clumped_snps)
        prep_pQTL_data(protname=args.protname,dataset=args.dataset,pQTL_path=args.pQTL_path,output_header=args.output_header,clumped_snps=args.clumped_snps)

    else:
        # If an unknown command is provided, show the help message
        parser.print_help()

if __name__ == "__main__":
    main()
