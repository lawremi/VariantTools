#!/bin/bash
source /gne/research/apps/modules/common/bashrc
module load apps/igis

zcat /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN//SAM587376/SAM587376_sample_specific_variants.vcf.gz |variant_effect_predictor.pl --sift=b --polyphen=b --condel=b --hgnc --check_existing --output_file /gnet/is3/research/data/bioinfo/ngs_analysis/degenhj2/CGP3.0_colon_var_rerun/TN/varEffOut/SAM587376_variant_annotations.txt