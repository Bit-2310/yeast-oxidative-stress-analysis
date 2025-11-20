# Biological Insights From the Yeast Oxidative Stress Dataset

The latest Nextflow run compares two controls (SRR636633–34) against four hydrogen-peroxide treated replicates (SRR636635–38). edgeR reports 783 genes with _p_ < 0.05 (679 up, 104 down), capturing the canonical Yap1/Skn7 oxidative-stress signature.

## Strongly Induced Antioxidant and Stress Genes

| Systematic ID | Common Name | log2FC | Biological role |
| --- | --- | --- | --- |
| YCR021C | **HSP30** | 7.15 | Proton-pump stabilizing heat shock protein that shields the plasma membrane during ROS insults. |
| YOL052C-A | **DDR2** | 7.01 | Yap1/Skn7 target induced by DNA damage and peroxide; supports stress-tolerant cell walls. |
| YPL171C | **OYE3** | 6.67 | FMN-dependent oxidoreductase that detoxifies lipid-derived peroxides. |
| YDR453C | **GTT2** | 7.57 | Glutathione transferase that conjugates electrophiles produced during oxidative bursts. |
| YKL086W | **SRX1** | 5.92 | Sulfiredoxin that reactivates hyperoxidized peroxiredoxins (Tsa1/Tsa2). |
| YBR244W | **GPX2** | 3.67 | Cytosolic glutathione peroxidase that removes H₂O₂ and lipid hydroperoxides. |

These transcripts highlight the detoxification toolkit: glutathione cycling, peroxiredoxin repair, and membrane fortification.

## Metabolic Rebalancing and Growth Suppression

| Systematic ID | Common Name | log2FC | Biological role |
| --- | --- | --- | --- |
| Q0158 | 15S mitochondrial rRNA | -2.66 | Lowered mitochondrial ribosome production reduces oxidative phosphorylation. |
| Q0020 | **COX3** | -2.67 | Cytochrome c oxidase subunit III; repression limits electron-transport chain flux. |
| LSR1 | U3 snoRNA | -1.50 | Required for 18S rRNA maturation, reflecting reduced ribosome biogenesis. |
| YFR032C | **RPL29** | -1.82 | 60S ribosomal protein; decreased abundance mirrors translational slowdown. |

The stress program therefore reallocates energy toward ROS mitigation while throttling mitochondrial respiration and ribosome assembly—an expected trade-off during acute peroxide exposure.
