library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(MRPRESSO)

args = commandArgs(trailingOnly=TRUE)

output_dir = 'protein_to_T2D/'

# read the outcome GWAS
protname <- args[1]
dataset = args[2]
N_PQTL = as.integer(args[3])
N_GWAS = as.double(args[4])

outcome_path = paste0(protname,'_',dataset,'_chrALL_clumped_withmeta_GWAS.txt')
# outcome_path <- '/scratch/richards/public/decode_proteomics_2021/4337_49_CRP_CRP.txt.gz'
# protname <- '4337_49_CRP_CRP'
outcome_GWAS <- vroom(outcome_path)

#exposure
exp_path = paste0(protname,'_',dataset,'_chrALL_clumped_withmeta.txt')
exp_dat <- read_exposure_data(filename = exp_path,
                              sep='\t',
                              snp_col='RSID',
                              beta_col = 'Effect',
                              se_col = 'StdErr',
                              effect_allele_col = 'EffectAllele',
                              other_allele_col = 'nonEffectAllele',
                              eaf_col = 'Freq1',
                              pval_col = 'P')

# read the outcome
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$RSID, type="outcome", snp_col="RSID", beta_col="Effect", se_col = "StdErr", eaf_col = "Freq1", 
                                 effect_allele_col = "EffectAllele", other_allele_col = "nonEffectAllele", pval_col = "P", chr_col = "CHROM", pos_col = "POS")
# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)
# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))

# mr result
mr_results <- mr(exp_dat_outcome)
mr_name <- paste0(output_dir, protname, ".",dataset,".mr.txt")
write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, protname, ".",dataset,".or.txt")
write.table(OR, file=OR_name, sep = '\t', quote = F, row.names = F)

# scatter plot
pdf_name <- paste0(output_dir, protname, ".",dataset,".scatter.pdf")
pdf(pdf_name, width = 10)
mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1) + theme_minimal()
dev.off() 

# horizontal pleiotropy
tryCatch({ 
  pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
  pleio_name <- paste0(output_dir,protname, ".",dataset,".pleio.txt")
  write.table(pleio_res, file=pleio_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# hetero test
tryCatch({ 
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$isquared <- abs(100*(hetero_res$Q - hetero_res$Q_df)/hetero_res$Q)  # I2 = 100%×(Q - df)/Q
  hetero_name <- paste0(output_dir,protname, ".",dataset,".hetero.txt")
  write.table(hetero_res, file=hetero_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# steiger
print('RUNNING STEIGER')
exp_dat_outcome$samplesize.exposure <- N_PQTL
exp_dat_outcome$samplesize.outcome <- N_GWAS
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, protname, ".",dataset,".steiger.txt")
write.table(steiger, file=steiger_name, sep = '\t', quote = F, row.names = F)

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

