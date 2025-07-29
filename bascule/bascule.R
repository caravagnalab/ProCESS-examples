options(bitmapType='cairo')
library(bascule)
library(optparse)
library(patchwork)
library(dplyr)

source("../getters/process_getters.R")
source("../getters/tumourevo_getters.R")
reticulate::use_condaenv("/orfeo/scratch/cdslab/ggandolfi/miniconda/envs/bascule-env")
py = reticulate::import("pybascule")


option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--purity"), type = "character", default = '0.9'),
                    make_option(c("--coverage"), type = "character", default = '50'),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn <- opt$spn_id #"SPN01"
coverage <- as.numeric(opt$coverage)
purity <- as.numeric(opt$purity)
cna_caller <- opt$cna_caller
vcf_caller <- opt$vcf_caller


data_dir = '/orfeo/scratch/cdslab/shared/SCOUT'

contexts <- c("SBS96","ID83")

reference_cat = list("SBS"=COSMIC_sbs_filt, "ID"=COSMIC_indels)
counts <- list()
spns <- list(spn)

counts_all <- list()
for (context in contexts){
  counts <- list()
  for (spn in spns){
    signature_class = gsub('[[:digit:]]+', '', context)
    
    sigprofiler_results <- get_tumourevo_signatures(spn = spn,coverage = coverage,purity = purity,
                                                    tool = "SigProfiler",context = context,
                                                    vcf_caller = vcf_caller,cna_caller = cna_caller)
    sig_matrix <- read.table(sigprofiler_results$context_matrix,header = T)
    rownames(sig_matrix) <- sig_matrix$MutationType
    sig_matrix <- sig_matrix %>% select(!MutationType)
    sig_counts <- t(sig_matrix) %>% as.data.frame()
    rownames(sig_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(sig_counts))
    counts[[spn]] <- sig_counts
  }
  counts_all[[context]] <- do.call("rbind",counts)
}
counts <- list("SBS"=counts_all$SBS96,"ID"=counts_all$ID83)
# reference_cat = list("SBS"=COSMIC_sbs_filt, "DBS"=COSMIC_dbs)


x = bascule::fit(
  counts = counts,
  k_list = 0:7,
  cluster = NULL,
  reference_cat = reference_cat,
  keep_sigs = c("SBS1","SBS5"),
  hyperparameters = NULL,
  lr = 0.005,
  optim_gamma = 0.1,
  n_steps = 3000,
  py = py,
  enumer = "parallel",
  nonparametric = TRUE,
  autoguide = FALSE,
  filter_dn = FALSE,
  min_exposure = 0.1,
  CUDA = TRUE,
  compile = FALSE,
  store_parameters = FALSE,
  store_fits = TRUE,
  seed_list = 10
)


x_refined = refine_denovo_signatures(x)
#x_refined_cluster = fit_clustering(x_refined, cluster=length(spns))
#x_refined_cluster = merge_clusters(x_refined_cluster)

output_dir <- file.path(data_dir,spn,"tumourevo",paste0(coverage,"x_",purity,"p_",vcf_caller,"_",cna_caller),"signature_deconvolution","BASCULE","SCOUT")
dir.create(path = output_dir,recursive = T)
saveRDS(object = x_refined,file = file.path(output_dir,"bascule_refined_fit.rds"))
saveRDS(object = x,file = file.path(output_dir,"bascule_fit.rds"))
################################################################################

exp = plot_exposures(x=x_refined)
sig = plot_signatures(x=x_refined)
plot_fit = exp/sig
ggplot2::ggsave(filename=file.path(output_dir,"bascule_fit.png"),plot=plot_fit)






