clc; clear all; close all; beep off;

% Parameters

p.gamma = 1/9;
p.alpha = 6;

% controls

N = 10;

scl1 = 0.1;
scl2 = 0.1;


%--------------------------------------------------------------------------
% Solve for P 
%--------------------------------------------------------------------------


mu1_u = sqrt(p.alpha)
mu2_u = 1/sqrt(p.gamma)

P = zeros(4,N+1,N+1);
P(:,1,1) = [1;0;0;0];
P(:,2,1) = scl1*[ 1; mu1_u; 0; 0];
P(:,1,2) = scl2*[0; 0; 1; mu2_u];

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

















%--------------------------------------------------------------------------
% Error testing of the parametrization method
%--------------------------------------------------------------------------

% decay rate of coefficients based on scaling of eigenvectors

figure 
hold on
coefAxis = 0:N;
plot(coefAxis, log(coeffMag)/log(10), 'b.','MarkerSize',18)


pnts = 11;
sigma1 = linspace(-1,1,pnts);
sigma2 = linspace(-1,1,pnts);

res_err = zeros(length(sigma1),length(sigma2));
for j = 1:length(sigma1)
    sig1 = sigma1(j);
    for k = 1:length(sigma2)
        sig2 = sigma2(k);
       
        y0 = eval_P(P,N,N,sig1,sig2);
       
        res_err(j,k) = norm(eval_P_test(P,N,N,sig1,sig2,mu1_u,mu2_u)- ...
            profile_ode(0,eval_P(P,N,N,sig1,sig2),p)); 
    end
end

max_res_err = max(max(res_err))


% return
%--------------------------------------------------------------------------
% Solve the profile
%--------------------------------------------------------------------------

% dependent controls
s.I = 10;
s.L = -s.I;
s.R = s.I;

% solve the profile
ode_fun = @(x,y,params)(profile_ode(x,y,p));
bc_fun = @(ya,yb,params)(profile_bc(ya,yb,p));
x = linspace(s.L,0,30);
guess = @(x)(profile_guess(x,p)); % NOT USING CONTINUATION
solinit = bvpinit(x,guess);
options = bvpset('RelTol',1e-6,'AbsTol',1e-6);
lastwarn('');
s.sol = bvp5c(ode_fun,bc_fun,solinit,options);
if strcmp('',lastwarn())==0
    error('warning given in solving the boundary value problem');
end

if length(s.sol.x) > 1000
   error('too many mesh points in profile solution'); 
end
dom = linspace(s.L,0,200);
end_state = deval(s.sol,s.L);
max_diff = 0;
for j = 1:length(dom)
    max_diff = max(max_diff,norm(end_state-deval(s.sol,dom(j))));
end
if max_diff < 1e-5
   error('solution is constant'); 
end

% % plot the profile
% plot_profile(p,s);

%--------------------------------------------------------------------------
% first order approximation
%--------------------------------------------------------------------------

L = -1.8;

z0 = deval(s.sol,L);

pnts = 11;
sigma1 = linspace(-1,1,pnts);
sigma2 = linspace(-1,1,pnts);
diff = zeros(length(sigma1),length(sigma2));
for j = 1:length(sigma1)
    for k = 1:length(sigma2)
        sig1 = sigma1(j);
        sig2 = sigma2(k);
        
        diff(j,k) = norm(eval_P(P,N,N,sig1,sig2)-z0);
                
    end
end

min_row = zeros(1,pnts);
for j = 1:length(min_row)
   min_row(j) = find(diff(j,:) == min(diff(j,:))); 
end

%--------------------------------------------------------------------------
% boundary value problem to find sigma1 and sigma2
%--------------------------------------------------------------------------


options = bvpset('AbsTol',10^(-8), 'RelTol',10^(-6));

guess = @(x)(interpolate(x,s.sol));
bc = @(ya,yb,var)(bc_parametrization(ya,yb,var,P,N));
x_grid = linspace(L,0,30);
var_guess = [0;0];
ode_bvp = @(x,y,var)(ode_paramtrization(x,y,var,p));
solinit = bvpinit(x_grid,guess,var_guess);
sol = bvp5c(ode_bvp,bc,solinit,options);
var = sol.parameters;

var

figure;
plot(sol.x,sol.y)

diff = eval_P(P,N,N,var(1),var(2))-deval(sol,L)


return
%%%%%%%%%%%%%%%%%%
figure 
hold on
coefAxis = 0:N;
plot(coefAxis, log(coeffMag)/log(10), 'b.','MarkerSize',18)

sigma1 = linspace(-0.8,0.8,11);
sigma2 = linspace(-0.8,0.8,11);


T = -0.0001;
err = zeros(length(sigma1),length(sigma2));
for j = 1:length(sigma1)
    for k = 1:length(sigma2)
        sig1 = sigma1(j);
        sig2 = sigma2(k);
        

        
        y0 = eval_P(P,N,N,sig1,sig2);
        
        ode = @(x,y)(profile_ode(x,y,p));
        options = odeset('RelTol',1e-13,'AbsTol',1e-12);
        sol = ode15s(ode,[0,T],y0,options);
        err(j,k) = norm(deval(sol,sol.x(end))- ...
        eval_P(P,N,N,sig1,sig2*exp(mu2_u*T)));
        
    end
end


max(err)

        ode = @(x,y)(profile_ode(x,y,p));


err2 = zeros(length(sigma1),length(sigma2));
for j = 1:length(sigma1)
    for k = 1:length(sigma2)
        sig1 = sigma1(j);
        sig2 = sigma2(k);
        

        
        y0 = eval_P(P,N,N,sig1,sig2);
        
        
        
        
        err2(j,k) = norm(eval_P_test(P,N,N,sig1,sig2,mu1_u,mu2_u)- ...
            ode(0,eval_P(P,N,N,sig1,sig2)));
        
    end
end






