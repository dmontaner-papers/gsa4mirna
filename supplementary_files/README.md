
This folder contains supplementary materials for the paper entitled "Integrated Gene Set Analysis for microRNA Studies" by Garcia-Garcia F, Panadero J, Dopazo J and Montaner D.


File content
------------

- res_edger_paired ### .xlsx   :  results from the __paired__   miRNA differential expression analysis using edgeR.
- res_edger_unpaired ### .xlsx :  results from the __unpaired__ miRNA differential expression analysis using edgeR.

- transfer_index_paired ### .xlsx   :  gene transferred statistics for the __paired__   analysis.
- transfer_index_unpaired ### .xlsx :  gene transferred statistics for the __unpaired__ analysis.

- res_gsa_paired ### bp.xls   :  results from the gene set analysis of the Gene Ontology Biological Processes in the __paired__   comparison.
- res_gsa_unpaired ### bp.xls :  results from the gene set analysis of the Gene Ontology Biological Processes in the __unpaired__ comparison.

- res_gsa_paired ### cc.xls   :  results from the gene set analysis of the Gene Ontology Cellular Component in the __paired__   comparison.
- res_gsa_unpaired ### cc.xls :  results from the gene set analysis of the Gene Ontology Cellular Component in the __unpaired__ comparison.

- res_gsa_paired ### mf.xls   :  results from the gene set analysis of the Gene Ontology Molecular Function in the __paired__   comparison.
- res_gsa_unpaired ### mf.xls :  results from the gene set analysis of the Gene Ontology Molecular Function in the __unpaired__ comparison.

Where ### indicates each cancer type.


- common_enrichment_paired.xlsx   :  indicates GO terms commonly enriched in several cancer types for the __paired__ comparison.
- common_enrichment_unpaired.xlsx :  indicates GO terms commonly enriched in several cancer types for the __unpaired__ comparison.

- allStudyInfo.xlsx  :  summary information of each of the 20 cancer types analyzed.


Plot interpretation
-------------------

- paired_explore_transfer ### .png: exploring the transferred index for miRNA and gene for each __paired__ study.
- paired_size_effect ### .png: these plots show the relationship between the functional block size and its significance, evaluating p-values and a combined indicator (-log p-value * sign (log OR)), in __paired__ studies.
- unpaired_size_effect ### .png: these set of plots show the relationship between the functional block size and 
its significance, evaluating p-values and a combined indicator (-log p-value * sign (log OR)), in __unpaired__ studies.

Where ### indicates each cancer type.

- paired_cor @@@.png :    plots of inhibition effect correlation, between all __paired__ studies.
- unpaired_cor @@@.png :  plots of inhibition effect correlation, between all __unpaired__ studies.
- paired_dist_of_cor@@@.png:    inhibition effect correlation distribution for all __paired__ studies.
- unpaired_dist_of_cor@@@.png:  inhibition effect correlation distribution for all __unpaired__ studies.

Where @@@ indicates each ontology: BP (Biological Process), CC (Cellular Component), MF (Molecular Function).

- paired_cor_rindex#.png:    ranking index correlation for all __paired__ studies.
- unpaired_cor_rindex#.png:  ranking index correlation for all __unpaired__ studies.

Where #:  0 (ranking index correlation miRNA), N (ranking index correlation genes), T (transferred ranking index).

- paired_inhibition_effect_correlation_across_ontologies:   evaluating inhibition effect correlation across ontologies for all __paired__ studies.
- unpaired_inhibition_effect_correlation_across_ontologies: evaluating inhibition effect correlation across ontologies for all __unpaired__ studies.

- paired_rindex_boxplot.png:    comparing distributions of several ranking index for all __paired__ studies.
- unpaired_rindex_boxplot.png:  comparing distributions of several ranking index for all __unpaired__ studies.

- paired_rindex_cor_vs_cor.png:    correlation between cancer types at gene transferred level plotted against correlation at miRNA level; __paired__ studies.
- unpaired_rindex_cor_vs_cor.png:  correlation between cancer types at gene transferred level plotted against correlation at miRNA level; __unpaired__ studies.


- inhibition_effect_cor_paired_vs_unpaired_@@@.png: displays correlation between the paired and the unpaired results at gene set level.
- inhibition_effect_paired_vs_unpaired_###.png:     displays, for each cancer, scatter plots of paired vs. unpaired results at gene set level.


- gene_vs_mirna_level_gsa_###.png: display the correlation between the GSA analysis carried out at miRNA level (Godard's paradigm) and at gene level after "transference".
  Each dot represents a GO term. X and Y values are derived from p-values and signs of the log odds ratios resulting from the mdgsa analysis
  (similar to equation 1 of the paper but at GO level instead of at miRNA level). 
