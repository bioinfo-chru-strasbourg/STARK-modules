# from PBMC_example$bulk_DE_cors # 59 sex-specific DEGs in bulk PBMC (up-regulated = female-biased)
# https://cran.r-project.org/web/packages/scMappR/vignettes/scMappR_Vignette.html

# Small list of sex related gene
#male.features <- c("Sry", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Usp9y", "Zfy1", "Zfy2", "Rps4y1")
#female.features <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Utx", "Rps4x", "Jarid1c", "Prdx4")

# Extended list of sex related gene
#female.features <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Rps4x", "Jarid1c", "Prdx4", "Med14","Zfx", "Ocrl", "Ddx39b", "Gpm6b", "Fthl17", "Esr1", "Foxl2", "Cyp19a1", "Wnt4", "Amhr2", "Eda2r")
#male.features <- c("Sry", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Usp9y", "Zfy1", "Zfy2", "Rps4y1", "Smcy", "Eif1ay", "Ssty1", "Ssty2", "Rbmy1a1", "Sox30", "Sox9", "Amh", "Insl3", "Hsd17b3", "Cyp17a1")

# male.features <- c("Abca13", "Adam22", "Arhgap31-as1", "Bcorp1", "Bpi", "Camp", "Cdy4p", "Ceacam19","Ceacam6", "Ceacam8", "Celsr2", "Crisp3", "Ctd-2308n23.2", "Ddx3y", "Defa3", "Defa4","Eif1ay", "Fbn1", "Gstt1", "Inhba", "Kalp", "Kdm5d", "Kif9-as1", "Lcn2", "Linc00278","Ltf", "Mmp8", "Mpo", "Prky", "Rn7skp282", "Rp11-424g14.1", "Rp6-99m1.2", "Rps4y1", "Slc1a2", "Tmsb4y", "Ttty14", "Ttty15", "Txlng2p", "Uty", "Usp9y", "Xist", "Ythdf3-as1", "Zfy", "Zfy-as1", "Znf839p1")
# female.features <- c("Ac084018.1", "Ighg3", "Ighv1-46", "Igkv1-27", "Igkv3-15", "Igkv3d-11", "Iglv1-40", "Iglv1-47", "Iglv3-1", "Iglv6-57", "Ift81", "Kdm5c", "Prkx", "Tsix", "Xist")

# curated
female.features.mouse <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Rps4x", "Prdx4", "Med14","Zfx", "Ocrl", "Ddx39b", "Gpm6b", "Esr1", "Foxl2", "Wnt4", "Amhr2", "Eda2r")
male.features.mouse <- c("Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Sox30", "Sox9", "Amh", "Insl3")

female_features.human <- c("IGKV1-27", "IGHG3", "IGHV1-46", "IGKV3-15", "IGKV3D-11", "IGLV3-1", "IGLV1-47", "IFT81", "IGLV1-40")
male_features.human <- c("XIST", "TSIX", "DDX3Y", "PRKY", "KDM5D", "RPS4Y1", "UTY", "USP9Y", "EIF1AY", "TXLNG2P", "ZFY", "TTTY15", "BCORP1", 
"LINC00278", "RP11-424G14.1", "TMSB4Y", "LTF", "RN7SKP282", "CEACAM8", "ZNF839P1", "ABCA13", "ZFY-AS1", "TTTY14", "DEFA4", "CEACAM6", "DEFA3", "BPI", "IGLV6-57", "RP6-99M1.2", 
 "KALP", "CAMP", "CRISP3", "KDM5C", "MPO", "PRKX", "KIF9-AS1", "ADAM22", "ARHGAP31-AS1", "CELSR2", "CDY4P", "AC084018.1", 
"LCN2", "CTD-2308N23.2", "INHBA", "FBN1", "SLC1A2", "GSTT1", "CEACAM19", "MMP8", "YTHDF3-AS1")

#Y-linked Genes (Male-biased)

#    Canonical Y-linked genes like Ddx3y, Eif2s3y, Uty, and Kdm5d are well-known for their roles in male sex determination and spermatogenesis. These genes are commonly referenced in mouse gene expression studies and are essential for male development.
#        Source: Mouse Genome Informatics (MGI) database, 10x Genomics reference datasets, and studies related to Y chromosome gene expression.
#        Example paper: Bellott, D. W., et al. (2014). "A human Y-chromosome gene is retained widely across 16 mammals." Nature.

#    Other Y-linked genes, such as Sry, Zfy1, and Zfy2, are involved in sex determination and testis development.
#        Source: The Y Chromosome Consortium and Mouse ENCODE Project.

#    Spermatogenesis-related genes like Rbmy, Tspy, and Ssty1/2 are involved in sperm development.
#        Source: Research articles on spermatogenesis in male mice and specific Y-chromosome gene clusters.
#        Example paper: Reynard, L. N., & Turner, J. M. A. (2009). "Interactions of meiotic proteins with the Y chromosome." Chromosome Research.

#X-linked Genes (Female-biased or X-Chromosome Genes)

#    X-inactivation-related genes like Xist, Tsix, and Rlim play crucial roles in dosage compensation between sexes by silencing one X chromosome in female cells.
#        Source: Research on X-chromosome inactivation in mammals.
#        Example paper: Lee, J. T. (2011). "Gracefully aging at 50, X-chromosome inactivation becomes a paradigm for RNA and chromatin control." Nature Reviews Molecular Cell Biology.

#    Dosage-sensitive genes like Kdm6a and Eif2s3x, which escape X-inactivation, are expressed in both male and female cells, but often at higher levels in females due to their escape from inactivation.
#        Source: Studies on X-escape genes and expression balance between sexes.
#        Example paper: Tukiainen, T., et al. (2017). "Landscape of X chromosome inactivation across human tissues." Nature.

#    X-linked immune-related genes such as Tlr7 and Foxp3, which are involved in immune system function, are known to have sex-biased expression in female-biased autoimmunity.
#        Source: Immune-related gene expression studies in mice.
#        Example paper: Souyris, M., et al. (2018). "Toll-like receptor 7-dependent autoimmunity in a mouse model of lupus." Science Immunology.

#    Other X-linked genes, including Ddx3x, Usp9x, and Fmr1, are critical for various cellular functions and are known to show differential expression between the sexes.
#        Source: 10x Genomics datasets, MGI database, and various papers on X-chromosome gene expression in mouse.
#        Example paper: Berletch, J. B., et al. (2011). "Escape from X inactivation varies in mouse tissues." PLoS Genetics.

# Expanded list of mouse sex-related genes
# male.features.mouse = c(
#    "Ddx3y", "Eif2s3y", "Uty", "Kdm5d", "Zfy2", "Zfy1", "Rbm31y", "Sry", 
#   "Jarid1d", "Eif2s3y", "Tspy", "Tmsb4y", "Tssk1", "Tssk2", "Ssty1", "Ssty2",
#    "Usp9y", "Rbm3y", "Smcy", "Ube1y", "Zfy", "Rbmy", "Dby", "Rps4y1", "Gm29650",
#    "Gm10352", "Gm19240", "Cnbp2")
  
  # X-linked genes (female-biased or X-chromosome genes)
#  female.features.mouse = c(
#    "Xist", "Tsix", "Kdm6a", "Eif2s3x", "Utx", "Ddx3x", 
#    "Jarid1c", "Kdm5c", "Rlim", "Stag2", "Btk", "Hdac6", "Mecp2", "Sytl2", "Timp1", 
#    "Xiap", "Atp7a", "Akap14", "G6pdx", "Ogt", "Timm17b", "Ccdc22", "Fhl1", "Tlr7", 
#    "Scml2", "Shroom4", "Armcx1", "Armcx2", "Armcx3", "Armcx4", "Armcx5", "Armcx6",
#    "Bcl6", "Eda", "Iqsec2", "Map3k15", "Prps1", "Rpl10", "Sdha", "Slc9a6", "Usp27x",
#    "Zfx", "Fmr1", "Gpr143", "Tfe3", "Amelx", "Foxp3", "Usp9x")

#male.features.human <- c(
#  "DDX3Y", "EIF2S3Y", "UTY", "KDM5D", "ZFY", "RBMY1A1", "SRY",
#  "USP9Y", "RPS4Y1", "RPS4Y2", "ZFY2", "TMSB4Y", "TSPY1", "PRY", "VCY", 
#  "SLC25A6", "BPY2", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "TTTY14", "HSFY1"
#)

#female.features.human <- c(
#  "XIST", "TSIX", "KDM6A", "EIF2S3X", "UTX", "DDX3X", 
#  "KDM5C", "RLIM", "STAG2", "BTK", "HDAC6", "MECP2", "SYTL2", "TIMP1", 
#  "XIAP", "ATP7A", "G6PD", "OGT", "TIMM17B", "CCDC22", "FHL1", "TLR7", 
#  "SCML2", "SHROOM4", "ARMCX1", "ARMCX2", "ARMCX3", "ARMCX4", "ARMCX5", "ARMCX6",
#  "BCL6", "EDA", "IQSEC2", "PRPS1", "RPL10", "SDHA", "SLC9A6", "USP27X",
#  "ZFX", "FMR1", "GPR143", "TFE3", "AMELX", "FOXP3", "USP9X"
#)
