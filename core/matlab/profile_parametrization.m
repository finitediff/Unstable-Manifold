function [P,max_res_err] = profile_parametrization(p,N)

    
scl_left = 0.1;
scl_right = 2;
[P,order_lower,coeffMag] = solve_P(p,N,scl_left);
[P,order_upper,coeffMag] = solve_P(p,N,scl_right);

if order_lower > -17
   error('desired scaling may not be enclosed'); 
elseif order_upper < -14
   error('desired scaling may not be enclosed'); 
end

order = 1;

while order > -14 || order < -17
    scl_mid = 0.5*(scl_left+scl_right);
    [P,order] = solve_P(p,N,scl_mid);
    if order > -14
        scl_right = scl_mid;
    elseif order < -17
        scl_left = scl_mid;
    end
end


% decay rate of coefficients based on scaling of eigenvectors

% [P,order,coeffMag] = solve_P(p,N,scl_mid);
% 
% figure 
% hold on
% coefAxis = 0:N;
% plot(coefAxis, log(coeffMag)/log(10), 'b.','MarkerSize',18)


%--------------------------------------------------------------------------
% Error testing of the parametrization method
%--------------------------------------------------------------------------

mu1_u = sqrt(p.alpha);
mu2_u = 1/sqrt(p.gamma);

pnts = 11;
sigma1 = linspace(-1,1,pnts);
sigma2 = linspace(-1,1,pnts);

res_err = zeros(length(sigma1),length(sigma2));
for j = 1:length(sigma1)
    sig1 = sigma1(j);
    for k = 1:length(sigma2)
        sig2 = sigma2(k);
              
        res_err(j,k) = norm(eval_P_test(P,N,N,sig1,sig2,mu1_u,mu2_u)- ...
            profile_ode(0,eval_P(P,N,N,sig1,sig2),p)); 
    end
end

max_res_err = max(max(res_err));

% -------------------------------------------------------------------------
% get parametrization
%--------------------------------------------------------------------------

function [P,order,coeffMag] = solve_P(p,N,scl)

    mu1_u = sqrt(p.alpha);
    mu2_u = 1/sqrt(p.gamma);

    P = zeros(4,N+1,N+1);
    P(:,1,1) = [1;0;0;0];
    P(:,2,1) = scl*[ 1; mu1_u; 0; 0];
    P(:,1,2) = scl*[0; 0; 1; mu2_u];

    J = zeros(4);
    J(1,2) = 1;
    J(2,1) = p.alpha;
    J(3,4) = 1;
    J(4,3) = 1/p.gamma;

    I4 = eye(4);

    coeffMag = zeros(1, N+1);
    coeffMag(1) = norm(P(:, 1, 1), inf);
    coeffMag(2) = max([norm(P(:, 2, 1), inf), norm(P(:, 1, 2), inf)]);


    for order = 2:N
        thisMax = 0;
        for suborder =  0:order
            n = order-suborder;
            m = suborder;
            temp = cauchy_prod3_smaller_terms(P(1,1:m+1,1:n+1),P(3,1:m+1,1:n+1),P(3,1:m+1,1:n+1));
            P(:,m+1,n+1) = ((m*mu1_u+n*mu2_u)*I4-J)\[0;temp;0;-temp/p.gamma];
            thisNorm = norm(P(:, m+1, n+1), inf);
            if thisNorm > thisMax
                thisMax = thisNorm;
            end
        end
        coeffMag(order+1) = thisMax;
    end


    
    order = log(coeffMag(end))/log(10);











