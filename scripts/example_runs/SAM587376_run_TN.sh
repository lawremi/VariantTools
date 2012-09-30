#!/bin/bash
source /gne/research/apps/modules/common/bashrc
module load apps/R

Rscript ~/svn/VariantTools/scripts/run_TN_wrapper.R /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/SAM587376/results/SAM587376.filtered_variants_granges.RData /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/SAM587377/results/SAM587377.filtered_variants_granges.RData /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/SAM587377/results/SAM587377.variants_granges.RData /gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/SAM587377/results/SAM587377.coverage.RData /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN/SAM587376 SAM587376 cgp_3.0_colon_TN_rerun
