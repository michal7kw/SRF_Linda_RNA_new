---
created: 2025-09-26T11:01
updated: 2025-09-26
---
*   `ES` (Enrichment Score)
*   `nes` (Normalized Enrichment Score)
    *   A positive `nes` means the pathway's genes are enriched at the top of ranked list (i.e., associated with up-regulated genes in cluster of interest).
    *   A negative `nes` means the pathway's genes are enriched at the bottom of list (i.e., associated with down-regulated genes).
*   `NOM p-val` (Nominal p-value): The uncorrected p-value.
*   `fdr` (False Discovery Rate): `NOM p-val` adjusted for multiple hypothesis testing. 
*   `FWER p-val` (Family-Wise Error Rate p-value): A more conservative corrected p-value than `fdr`.
*   `Tag %`: The percentage of genes from the gene set that are in the "leading edge" subset.
*   `Gene %`: The percentage of genes in entire ranked list that belong to this specific gene set.
*   `Lead_genes` (Leading Edge Genes): The specific genes from the pathway that contributed most to the enrichment score.
*   `Gene_Set`: The name of the source database.
*   `Direction`: It is based on the `nes`:
    *   `Up`: Corresponds to a positive `nes`.
    *   `Down`: Corresponds to a negative `nes`.
*   `abs_nes` (Absolute `nes`)