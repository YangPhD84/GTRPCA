# ITRPCA 
To mine latent association features in multiple similarities and association data, we present an improved tensor robust principal component analysis (ITRPCA) method. First, we integrate the prior information of drug and disease to compute five indicators for drug similarity and two indicators for disease similarity. Considering that validated drug-disease associations are extremely sparse, a weighted k-nearest neighbor (WKNN) preprocessing step is employed to enrich the association matrix that aids in prediction. Then, we construct a drug tensor and a disease tensor by using multi-similarity matrices and an updated association matrix. Finally, we apply ITRPCA to isolate the low-rank tensor and noise information in these two tensors, respectively. We focus on the drug-disease association pairs in the clean low-rank tensor to infer promising indications for drugs.

# Requirements
* Matlab >= 2021

# Installation
ITRPCA can be downloaded by
```
git clone https://github.com/YangPhD84/ITRPCA
```
Installation has been tested in a Windows platform.

# Dataset Description
* DrugSimMat.mat: This file contains five drug similarity matrices;
* DiseaseSimMat.mat: This file contains two disease similarity matrices;
* Didr.mat: the disease-drug association's matrix;
* DrugsName: the DrugBank IDs of drugs;
* DiseasesName: the OMIM IDs of diseases;
* drug_ChemS: chemical structure similarity matrix;
* drug_AtcS: drug's ATC code similarity matrix;
* drug_SideS: side-effect similarity matrix;
* drug_DDIS: drug-drug interaction similarity matrix;
* drug_TargetS: drug's target profile similarity matrix;
* disease_PhS: disease phenotype similarity matrix;
* disease_DoS: disease ontology similarity matrix.

# Functions Description
* WKNN.m: this function can implement the WKNN preprocessing algorithm;
* fITRPCA.m: this function is the core algorithm of the ITRPCA model for implementing correlation prediction;
* ffindw.m: this function can implement the solution of the weight vector omega;
* itrpca_tnn_lp_stop.m: this function solves the problem of tensor robust principal component analysis based on the weighted tensor Schatten p-Norm;
* prox_l1.m: this function is the M-ADMM algorithm;
* prox_tnn.m: this function is used to update the low-rank tensor in the model;
* solve_Lp_w.m: this function can solve our ITRPCA model;
* Tprod.m: this function is capable of extracting the correlation matrix within a low-rank tensor.

# Instructions
We provide detailed step-by-step instructions for running ITRPCA model.
**Step 1**: add datasets\functions paths
```
addpath('Gold standard dataset');
addpath('Functions');
```
**Step 2**: load datasets with association matirx and similarity matrices

```
load DrugSimMat
load DiseaseSimMat
load didr
A_DR = didr;
[dn,dr] = size(A_DR);
Trr = zeros(dr,dr,1);
Trr(:,:,1) = drug_ChemS;
Trr(:,:,2) = drug_AtcS;
Trr(:,:,3) = drug_SideS;
Trr(:,:,4) = drug_DDIS;
Trr(:,:,5) = drug_TargetS;
Tdd = zeros(dn,dn,1);
Tdd(:,:,1) = disease_PhS;
Tdd(:,:,2) = disease_DoS;
```

# Installation
HGIMC can be downloaded by
```
git clone https://github.com/BioinformaticsCSU/HGIMC
```

We propose a new scheme with improved tensor robust principal component analysis (ITRPCA) to predict promising drug-associated indications. The code in this package implements ITRPCA for computational drug repositioning, which is implemented in Matlab.

# Contact:  
If you have any questions or suggestions with the code, please let us know. Contact Mengyun Yang at mengyun_yang@126.com
