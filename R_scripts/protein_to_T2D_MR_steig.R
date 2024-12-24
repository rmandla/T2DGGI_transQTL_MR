library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(MRPRESSO)

args = commandArgs(trailingOnly=TRUE)

output_dir = 'protein_to_T2D/'

# read the outcome GWAS
#protname <- args[1]
dataset = args[1]
N_PQTL = as.integer(args[2])
N_GWAS = as.double(args[3])
protname = 'IGFBP2'

outcome_path = paste0('IGFBP2_',dataset,'_GWAS_harmonized.txt')
exp_path = paste0('IGFBP2_',dataset,'_harmonized.txt')
# outcome_path <- '/scratch/richards/public/decode_proteomics_2021/4337_49_CRP_CRP.txt.gz'
# protname <- '4337_49_CRP_CRP'
outcome_GWAS <- vroom(outcome_path)

#exposure
#exp_path = paste0(protname,'_',dataset,'_chrALL_clumped_withmeta.txt')
exp_dat <- read_exposure_data(filename = exp_path,
                              sep='\t',
                              snp_col='rsid',
                              beta_col = 'Effect',
                              se_col = 'SE',
                              effect_allele_col = 'eAllele',
                              other_allele_col = 'oAllele',
                              eaf_col = 'AF',
                              pval_col = 'P')

# read the outcome
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$rsid, type="outcome", snp_col="rsid", beta_col="Effect", se_col = "SE", eaf_col = "AF", 
                                 effect_allele_col = "eAllele", other_allele_col = "oAllele", pval_col = "P", chr_col = "chr", pos_col = "pos")
# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)

exp_dat_outcome$samplesize.outcome <- N_GWAS
exp_dat_outcome$samplesize.exposure <- N_PQTL

exp_dat_outcome <- steiger_filtering(exp_dat_outcome)
exp_dat_outcome <- subset(exp_dat_outcome, steiger_dir & steiger_pval<0.05)

#if(nrow(exp_dat_outcome)>0){
# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))

# mr result
mr_results <- mr(exp_dat_outcome)
mr_name <- paste0(output_dir, protname, ".",dataset,".steig.mr.txt")
write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, protname, ".",dataset,".steig.or.txt")
write.table(OR, file=OR_name, sep = '\t', quote = F, row.names = F)

# scatter plot
pdf_name <- paste0(output_dir, protname, ".",dataset,".steig.scatter.pdf")
pdf(pdf_name, width = 10)
mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1) + theme_minimal()
dev.off() 

# horizontal pleiotropy
tryCatch({ 
  pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
  pleio_name <- paste0(output_dir,protname, ".",dataset,".steig.pleio.txt")
  write.table(pleio_res, file=pleio_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# hetero test
tryCatch({ 
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$isquared <- abs(100*(hetero_res$Q - hetero_res$Q_df)/hetero_res$Q)  # I2 = 100%×(Q - df)/Q
  hetero_name <- paste0(output_dir,protname, ".",dataset,".steig.hetero.txt")
  write.table(hetero_res, file=hetero_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# steiger
print('RUNNING STEIGER')
exp_dat_outcome$samplesize.exposure <- N_PQTL
exp_dat_outcome$samplesize.outcome <- N_GWAS
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, protname, ".",dataset,".steig.steiger.txt")
write.table(steiger, file=steiger_name, sep = '\t', quote = F, row.names = F)
#}
# # MR-PRESSO
# tryCatch({ 
#   mrpresso_res <- mr_presso(BetaExposure = "beta.exposure",
#                             BetaOutcome = "beta.outcome",
#                             SdOutcome = "se.outcome",
#                             SdExposure = "se.exposure",
#                             OUTLIERtest = TRUE,
#                             DISTORTIONtest = TRUE,
#                             data = exp_dat_outcome,
#                             NbDistribution = 15000,
#                             SignifThreshold = 0.05)
#   saveRDS(mrpresso_res, file = 'output_COL6A3/mrpresso.RDS')
#   
#   mrpresso_name <- paste0(output_dir, "mrpresso/", protname, ".mrpresso.txt")
#   
#   # add global_rss, global_pval, distortion_indices, distortion_coef, distortion_pval
#   mrpresso_df <- as.data.frame(mrpresso_res$`Main MR results`)
#   mrpresso_df %<>% mutate(#global_rss = mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs,
#     global_pval =  mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue,
#     #distorition_indices = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,
#     #distortion_coef = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
#     distortion_pval = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue)
#   
#   write_tsv(mrpresso_df, file =  mrpresso_name)      
#   
#   # additional information
#   mrpresso_global_rss <- mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs
#   mrpresso_global_pval <- mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue
#   mrpresso_distortion_indices <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
#   mrpresso_distortion_coef <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
#   mrpresso_distortion_pval <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
#   
#   mrpresso_add_res <- data.frame(mrpresso_global_rss = mrpresso_global_rss,
#                                  mrpresso_global_pval = mrpresso_global_pval,
#                                  mrpresso_distortion_indices = mrpresso_distortion_indices,
#                                  mrpresso_distortion_coef = mrpresso_distortion_coef,
#                                  mrpresso_distortion_pval = mrpresso_distortion_pval
#   )
#   
#   mrpresso_add_name <- paste0(output_dir, "mrpresso_add/", protname, ".mrpresso_add.txt")
#   write.table(mrpresso_add_res, file=mrpresso_add_name, sep = '\t', quote = F, row.names = F)
#   
# }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

