function [L,S,iter,stop1,stop2] = itrpca_tnn_lp_stop(X, lambda, weight, p, dimdim)
% Solve the Tensor Robust Principal Component Analysis based on Weighted Tensor Schatten p-Norm problem by ADMM
% min_{L,S} ||L||_sp^p+lambda*||S||_1, s.t. X=L+S

rho = 1.1;
mu = 1e-2;
max_mu = 1e10;
maxiter=5;
tol1=1e-3;
tol2=1e-4;

dim = size(X);
L = zeros(dim);
S = L;
Y = L;

stop1 = 1;
X_0=X;
for iter = 1 : maxiter
    Lk = L;
    Sk = S;
    % Update L
    [L,tnnL] = prox_tnn(-S+X-Y/mu,weight/mu,p);
    L(L<0)=0; % A bounded constrain
    L(L>255)=255;
    
    % Update S
    S = prox_l1(-L+X-Y/mu,lambda/mu);
    
    dY = L+S-X;
    
    % Coverge condition
    X_1=L;
    stop1_0 = stop1;
    sum_norm=0;
    for j=1:dimdim
        sum_norm=sum_norm+norm(X_1(:,:,j)-X_0(:,:,j),'fro')/norm(X_0(:,:,j),'fro');
    end
    stop1=sum_norm;
    stop2=abs(stop1-stop1_0)/max(1,abs(stop1_0));
    X_0=X_1;

    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);
    
    if (stop1 < tol1) && (stop2 < tol2)
        break
    end
end

end