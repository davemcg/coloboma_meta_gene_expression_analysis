# coloboma_meta_gene_expression_analysis

# Script run order
1. [src/microarray_workup.R]() 
  - creates [data/microarray_table.Rdata]() 
2. [src/NGS_processing.R]()
  - creates [data/NGS_processing.Rdata]()
3. [src/merge_microarray_NGS_normalize.R]()
  - creates [data/microarray_NGS_objects.Rdata]()
4. [src/src/limma_diff_testing.R]()
  - makes top tables and sva_counts
  
# Run

Rscript src/microarray_workup.R
Rscript src/NGS_processing.R
Rscript src/merge_microarray_NGS_normalize.R
Rscript src/limma_diff_testing.R
