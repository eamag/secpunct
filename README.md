Work on this research: https://docs.google.com/document/d/15fnkS3t_-PSZ33vVn_vBZBAA3nVtAvfqIvDHw5S6GVU/edit
### TODO:
- [x] Read databases in pandas, look up for eda in bioinf
- [x] Lift over hg18 -> hg19
- [x] 10-70bp (because of sites for rna polymerase) for make10bp in preprocessing.py
- [x] Cluster sec structs with true labels vs non-relevant
- [x] Binary LightGBM to predict relevant sec.struct (different phys. properties)
- [x] Repeat all for quadruplexes
- [ ] Repeat all for different chromosomes
- [x] Check logreg for feature importance and r-score instead of roc-auc