function [V1,V2,err1,err2] = evans_parametrization(p,N,P,lambda)

%--------------------------------------------------------------------------
% Solve for Evans function basis
%--------------------------------------------------------------------------

mu1_lam = sqrt(lambda+p.alpha);
mu2_lam = sqrt(lambda+1/p.gamma);

xi1_lam = [1;mu1_lam;0;0];
xi2_lam = [0;0;1;mu2_lam];

[V1,err1] = get_v(mu1_lam,xi1_lam,p,N,P,lambda);
[V2,err2] = get_v(mu2_lam,xi2_lam,p,N,P,lambda);

%--------------------------------------------------------------------------
% get parametrization v
%--------------------------------------------------------------------------


function [VP,max_res_err] = get_v(mu_lam,xi_lam,p,N,P,lambda)


mu1_u = sqrt(p.alpha);
mu2_u = 1/sqrt(p.gamma);

I4 = eye(4);

%
% Create the Bmn matrices
% 

Bmn = zeros(4,4,N,N);
Bmn(1,1,1,1) = -mu_lam;
Bmn(1,2,1,1) = 1;
Bmn(2,1,1,1) = lambda+P(3,1,1)^2+p.alpha;
Bmn(2,2,1,1) = -mu_lam;
Bmn(2,3,1,1) = 2*P(1,1,1)*P(3,1,1);
Bmn(3,3,1,1) = -mu_lam;
Bmn(3,4,1,1) = 1;
Bmn(4,1,1,1) = -P(3,1,1)^2/p.gamma;
Bmn(4,3,1,1) = lambda+(1-2*P(1,1,1)*P(3,1,1))/p.gamma;
Bmn(4,4,1,1) = -mu_lam;

for order = 1:N
    for m =  0:order
        n = order-m;
        
        temp1 = cauch_prod_dim2(squeeze(P(3,1:m+1,1:n+1)), ...
            squeeze(P(3,1:m+1,1:n+1)));
        
        temp2 = 2*cauch_prod_dim2(squeeze(P(1,1:m+1,1:n+1)), ...
            squeeze(P(3,1:m+1,1:n+1)));
        
        Bmn(2,1,m+1,n+1) = temp1;
        
        Bmn(4,1,m+1,n+1) = -temp1/p.gamma;
        
        Bmn(2,3,m+1,n+1) = temp2;
        
        Bmn(4,3,m+1,n+1) = -temp2/p.gamma;
        
    end
end


VP = zeros(4,N+1,N+1);
VP(:,1,1) = xi_lam;

B00 = squeeze(Bmn(:,:,1,1));

for order = 1:N
    for suborder =  0:order
        n = order-suborder;
        m = suborder;
        
        delta = ones(m+1,n+1);
        delta(m+1,n+1) = 0;
        summation = zeros(4,1); 
        for j = 0:m
            for k = 0:n
                summation = summation + delta(j+1,k+1)*squeeze(Bmn(:,:,m-j+1,n-k+1))* ...
                    squeeze(VP(:,j+1,k+1));
            end
        end
        
        VP(:,m+1,n+1) = ((m*mu1_u+n*mu2_u)*I4-B00)\summation;

        
    end
end

%
% Check the correctness of Evan manifold solution
% 

pnts = 11;
sigma1 = linspace(-1,1,pnts);
sigma2 = linspace(-1,1,pnts);

res_err = zeros(length(sigma1),length(sigma2));
for j = 1:length(sigma1)
    sig1 = sigma1(j);
    for k = 1:length(sigma2)
        sig2 = sigma2(k);
       
        u0 = eval_P(P,N,N,sig1,sig2);
        
        u = u0(1);
        v = u0(3);
        
        A0 = [0 1 0 0; 
             lambda+v^2+p.alpha 0 2*u*v 0;
             0 0 0 1;
             -v^2/p.gamma 0 lambda+(1-2*u*v)/p.gamma 0];
         
         y0 = eval_P(VP,N,N,sig1,sig2);
      
        res_err(j,k) = norm(eval_P_test(VP,N,N,sig1,sig2,mu1_u,mu2_u)- ...
            (A0 - mu_lam*I4)*y0);
    end
end

max_res_err = max(max(res_err));
















