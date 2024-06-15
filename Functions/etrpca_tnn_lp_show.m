function [L,S,obj,err,iter] = etrpca_tnn_lp_show(X, P_TMat, TestIds, lambda, weight, p, dimdim, opts)

% Solve the Tensor Robust Principal Component Analysis based on Weighted Tensor Schatten p-Norm problem by ADMM
%
% min_{L,S} ||L||_sp^p+lambda*||S||_1, s.t. X=L+S
%
% ---------------------------------------------
% Input:
%       X       -    d1*d2*d3 tensor
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       L       -    d1*d2*d3 tensor
%       S       -    d1*d2*d3 tensor
%       obj     -    objective function value
%       err     -    residual
%       iter    -    number of iterations
%
%
%
% Written by Pu Zhang

tol = 1e-18;
% tol = 1e-3;
max_iter = 50;%50
rho = 1.1;
mu = 1e-2;
max_mu = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBmUG = opts.DEBUG;          end

dim = size(X);
L = zeros(dim);
S = L;
Y = L;

for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % Update L
    [L,tnnL] = prox_tnn(-S+X-Y/mu,weight/mu,p);
    
    %%  每一步迭代都计算AUC值
             M_ResultMat123=L/255;
             [dn,dr] = size(P_TMat);
        for iii=1:dimdim
        M_ResultMat=M_ResultMat123((dim-dn+1):dim,1:dr,iii);
        [AUC_alone,Precision_alone,AUPR_alone]=Fun_Auc3(M_ResultMat,P_TMat,TestIds);
        fprintf('AUC=%4.4f,Precision=%4.4f,AUPR=%4.4f.\n',AUC_alone,Precision_alone,AUPR_alone)
        end
    %%
    % Update S
    S = prox_l1(-L+X-Y/mu,lambda/mu);
    
    dY = L+S-X;
    % Coverge condition
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    chg
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnL+lambda*norm(S(:),1);
            err = norm(dY(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                ', obj=' num2str(obj) ', err=' num2str(err)]);
        end
    end
    
    if chg < tol
        break;
    end
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);
    
    disp(['The ' num2str(iter)  '-th Iteration']);
end
obj = tnnL+lambda*norm(S(:),1);
err = norm(dY(:));