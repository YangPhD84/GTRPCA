clear all
addpath('Datasets');
addpath('Functions');

%% 1. Load Datesets
load dataset1_y2020
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

%% 2. GTRPCA algorithm
p = 0.9; 
K = 30;
rat1 = 0.1;
rat2 = 0.2;
A_recovery =fGTRPCA(Trr,Tdd,A_DR,p,K,rat1,rat2);




