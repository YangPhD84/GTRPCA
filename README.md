# ITRPCA 
To mine latent association features in multiple similarities and association data, we present an improved tensor robust principal component analysis (ITRPCA) method. First, we integrate the prior information of drug and disease to compute five indicators for drug similarity and two indicators for disease similarity. Considering that validated drug-disease associations are extremely sparse, a weighted k-nearest neighbor (WKNN) preprocessing step is employed to enrich the association matrix that aids in prediction. Then, we construct a drug tensor and a disease tensor by using multi-similarity matrices and an updated association matrix. Finally, we apply ITRPCA to isolate the low-rank tensor and noise information in these two tensors, respectively. We focus on the drug-disease association pairs in the clean low-rank tensor to infer promising indications for drugs.

* kdfvdj
```
dsafd
```
# Requirements
* Matlab >= 2014

# Installation
HGIMC can be downloaded by
```
git clone https://github.com/BioinformaticsCSU/HGIMC
```

We propose a new scheme with improved tensor robust principal component analysis (ITRPCA) to predict promising drug-associated indications. The code in this package implements ITRPCA for computational drug repositioning, which is implemented in Matlab.

# Contact:  
If you have any questions or suggestions with the code, please let us know. Contact Mengyun Yang at mengyun_yang@126.com
