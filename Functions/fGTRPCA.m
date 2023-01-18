function  M_ResultMat=fGTRPCA (Trr,Tdd,P_TMat,p,K,rat1,rat2)
%% WKNKN step
[dn,dr] = size(P_TMat);
P_TMat_new=WKNKN(P_TMat,mean(Tdd,3),mean(Trr,3),K,0.95);

%% Drug-tensor
[tr1,tr2,dr_num]=size(Trr);
Tdr=ones(dn,dr,dr_num);
for l1=1:dr_num
    Tdr(:,:,l1)=P_TMat_new;
end
R_ori=[Trr;Tdr]*255;

[n1,n2,n3] = size(R_ori);
n=min(n1,n2);
w = [];
for e1=1:dr_num
    ind1(e1) =ffindw(R_ori(:,:,e1), rat1);
end
a1=round(mean(ind1));
for e2=1:dr_num
    ind2(e2) =ffindw(R_ori(:,:,e2), rat2);
end
a2=-a1+2+round(mean(ind2));
w = [w; 1*ones(a1,1)];
w = [w; 2*ones(a2,1)];
w = [w; 4*ones(n-a1-a2,1)];

kao=1/(5*sqrt(n1*n2*n3));
[R_ResultMat123, E, iter,stop1,stop2]  = gtrpca_tnn_lp_stop(R_ori, kao, w, p,dr_num);
R_ResultMat123=R_ResultMat123/255;
R_ResultMat=mean(R_ResultMat123((n1-dn+1):n1,1:dr,1:dr_num),3);

%% Disease-tensor
[td1,td2,dd_num]=size(Tdd);
Tdr=ones(dn,dr,dd_num);
for l1=1:dd_num
    Tdr(:,:,l1)=P_TMat_new;
end
D_ori=[Tdr,Tdd]*255;

[nn1,nn2,nn3] = size(D_ori);
nn=min(nn1,nn2);
w = [];
for e1=1:dd_num
    indd1(e1) =ffindw(D_ori(:,:,e1), rat1);
end
b1=round(mean(indd1));
for e2=1:dd_num
    indd2(e2) =ffindw(D_ori(:,:,e2), rat2);
end
b2=-b1+2+round(mean(indd2));
w = [w; 1*ones(b1,1)];
w = [w; 2*ones(b2,1)];
w = [w; 4*ones(nn-b1-b2,1)];

kao=1/(5*sqrt(nn1*nn2*nn3));
[D_ResultMat123, E,iter,stop1,stop2]  = gtrpca_tnn_lp_stop(D_ori, kao, w, p,dd_num);
D_ResultMat123=D_ResultMat123/255;
D_ResultMat=mean(D_ResultMat123(1:dn,1:dr,1:dd_num),3);

%% Drug-disease association matrix
M_ResultMat=(R_ResultMat+D_ResultMat)/2;

