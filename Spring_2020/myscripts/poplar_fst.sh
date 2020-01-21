#!/bin/bash
vcftools --vcf /data/example_data/snp/poplar_437ind_51pop_129Ksnp.vcf \
--weir-fst-pop core.pop \
--weir-fst-pop edge.pop \
--out balsam_core_edge_fst
