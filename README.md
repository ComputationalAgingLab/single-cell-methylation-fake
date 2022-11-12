# single-cell-methylation-fake

This repository provides a set of evidences that the original procedure proposed in [Trapp et al.](https://www.nature.com/articles/s43587-021-00134-3) for inference epigenetic age in single cell methylation data contains mistakes and generates wrong results. I specifically considered one of the most strange results which is called "ground zero", i.e. achieving the global minimum in epigenetic age during embryogenesis. I show that the ground zero is actually an artifact of the inherently wrong procedure of maximum likelihood estimation (MLE) proposed by the authors. In short, authors proposed to find polynomial likelihood function by applying brute force algorithm manually defined interval of search while the global solution can be found. I show the global solution for MLE problem provides other results than the authors proposed in their paper. The concrete procedure is provided in the notebook `proof`.

# Acknowledgements

# Cite

