#!/bin/bash
source /gne/research/apps/modules/common/bashrc
module load apps/R

Rscript ~/svn/VariantTools/scripts/run_callVar_wrapper.R /gnet/is3/research/data/bioinfo/ngs_analysis/CGP_3.0/merged/NGS121/SAM587376/bams/SAM587376.analyzed.bam hg19_IGIS21 /gne/research/apps/ngs_pipeline/dev/x86_64-linux-2.6-sles11/share 75 /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/Exome/SAM587376 SAM587376 cgp_rerun
