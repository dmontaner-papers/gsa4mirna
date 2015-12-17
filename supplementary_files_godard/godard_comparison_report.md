% Integrated Gene Set Analysis for microRNA Studies
% Comparison of methods for Gene Set Analysis
% 2015-10-22



Methods
===============================================================

A GSA extension of Godard’s approach have been  computed straight forward using logistic regression
models and our Bioconductor library. This strategy will certainly retain Godard’s methodology good
characteristics while incorporating the benefits of the GSA approach over the ORA one.

We did find that the functional results at miRNA level (Godard’s generalization) and the ones at gene level (after transference as originally proposed)
have a significantly positive correlation. This indicates that overall both methodologies should provide
similar findings. The correlation is not very strong though as the methodologies are explicitly different.

Detailed results and scripts are available in [https://github.com/dmontaner-papers/gsa4mirna](https://github.com/dmontaner-papers/gsa4mirna)
(folders: supplementary_files_godard and scripts_godard).



Results
===============================================================

For each cancer type, there are several plots displaying the correlation between the GSA analysis carried out at miRNA level (Godard’s paradigm)
 and at gene level after “transference”. Each dot represents a GO term. 
 X and Y values are derived from p-values and signs of the log odds ratios resulting from the mdgsa analysis
  (similar to equation 1 of the paper but at GO level instead of at miRNA level).
\newpage

![BLCA](plots/gene_vs_mirna_level_gsa_blca.png)

![ESCA](plots/gene_vs_mirna_level_gsa_esca.png)

![KIRP](plots/gene_vs_mirna_level_gsa_kirp.png)

![PAAD](plots/gene_vs_mirna_level_gsa_paad.png)

![SKCM](plots/gene_vs_mirna_level_gsa_skcm.png)

![BRCA](plots/gene_vs_mirna_level_gsa_brca.png)

![HNSC](plots/gene_vs_mirna_level_gsa_hnsc.png)

![LIHC](plots/gene_vs_mirna_level_gsa_lihc.png)

![PCPG](plots/gene_vs_mirna_level_gsa_pcpg.png)

![STAD](plots/gene_vs_mirna_level_gsa_stad.png)

![CESC](plots/gene_vs_mirna_level_gsa_cesc.png)

![KICH](plots/gene_vs_mirna_level_gsa_kich.png)

![LUAD](plots/gene_vs_mirna_level_gsa_luad.png)

![PRAD](plots/gene_vs_mirna_level_gsa_prad.png)

![THCA](plots/gene_vs_mirna_level_gsa_thca.png)

![COAD](plots/gene_vs_mirna_level_gsa_coad.png)

![KIRC](plots/gene_vs_mirna_level_gsa_kirc.png)

![LUSC](plots/gene_vs_mirna_level_gsa_lusc.png)

![USEC](plots/gene_vs_mirna_level_gsa_ucec.png)

![READ](plots/gene_vs_mirna_level_gsa_read.png)
