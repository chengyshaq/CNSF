function [model_CNSF] = CNSF(X, Y, optmParameter)

 alpha          = optmParameter.alpha;

 gamma          = optmParameter.gamma;  
maxIter         = optmParameter.maxIter;
miniLossMargin  =  optmParameter.minimumLossMargin;
theta           = optmParameter.theta;

CL = Causal(Y') ;              
P = ICNE(Y);  
Y =P;
[~,num_dim]=size(X);
C = pdist2(Y', Y', 'cosine');
R = diag(sum(C)) - C;
R = CL.*R;

XTX = X'*X;
 W_s   = (XTX + theta*eye(num_dim)) \ (X'*Y);
 W_s_1 = W_s;
    iter    = 1;
    bk = 1;
    bk_1 = 1; 
    oldloss=0;
while iter <= maxIter
       Lip = sqrt(2*(norm(XTX)^2) +2*(2*alpha*norm(R)^2));
       W_s_k  = W_s + (bk_1 - 1)/bk * (W_s - W_s_1);
       Gw_s_k = W_s_k - 1/Lip * (-X'*Y+XTX *W_s_k+2*alpha*W_s_k*R);
       bk_1   = bk;
       bk     = (1 + sqrt(4*bk^2 + 1))/2;
       W_s_1 = W_s;
       W_s  = softthres(Gw_s_k,alpha/Lip);
        
        partA = 0.5*norm((X*W_s-Y),'fro')^2; 
        partB = trace(R*(W_s'*W_s));
        sparsity = sum(sum(W_s~=0));
        totalloss = partA + alpha*partB +  gamma*sparsity;

        if abs(oldloss - totalloss) <= miniLossMargin

             break;
        elseif totalloss <=0
             break;
        else
             oldloss = totalloss;
        end      
       iter=iter+1; 
end
    model_CNSF=W_s;
end
function W = softthres(W_t,gamma)
    W = max(W_t-gamma,0) - max(-W_t-gamma,0); 
end
