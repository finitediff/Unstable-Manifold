clc; clear all; close all; beep off;

%
% Parameters
%

p.gamma = 1/9;
p.alpha = 6;

%
% controls
%

N = 20; % number of terms in power expansion 

%--------------------------------------------------------------------------
% Solve for parametrization P of profile
%--------------------------------------------------------------------------
[P,max_res_err_of_profile] = profile_parametrization(p,N);

fprintf('\n\nMax residual error of parameterization of profile: %4.4g\n\n', ...
    max_res_err_of_profile);







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
% return

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

sig1 = var(1);
sig2 = var(2);

% figure;
% plot(sol.x,sol.y)
% 
% diff = eval_P(P,N,N,var(1),var(2))-deval(sol,L)

%--------------------------------------------------------------------------
% contour
%--------------------------------------------------------------------------


[s,e,m,c] = emcset(s,'front', LdimRdim(@A,s,p),'reg_reg_polar');
c.evans = @local_evans; % even Evans function
c.basisL = @analytic_basis_local;
c.basisR = @analytic_basis_local;


% 
% process and display data
%

% display a waitbar
c.stats = 'off';
c.refine = 'off';
c.ksteps = 0;
c.lambda_steps = 0;
c.tol = 0.2;
% m.options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Refine',1,'Stats','off','MaxOrder',2);
% m.ode_fun = @ode15s;

R = 1000; 
circpnts=30; imagpnts=30; spread=2; inner_pnts = 10; inner_radius = 0.01;
preimage=semicirc2(circpnts,imagpnts,inner_pnts,c.ksteps,R, ...
    spread,inner_radius,c.lambda_steps);

%--------------------------------------------------------------------------
% Evans function with parametrization and bvp solver
%--------------------------------------------------------------------------


E = zeros(size(preimage));

t1 = tic;
parfor j=1:length(E)
    E(j) = evans_parametrized(s,p,N,P,preimage(j),sig1,sig2,L);
end
time1 = toc(t1)

E = E/E(1);
w = [E,fliplr(conj(E))];
plot_evans(w)
winding_number(w)

E1 = E;

%--------------------------------------------------------------------------
% Evans function with bvp solver
%--------------------------------------------------------------------------

E = zeros(size(preimage));


t1 = tic;
parfor j=1:length(E)
    E(j) = evans_bvp(s,p,N,P,preimage(j),sig1,sig2,s.L);
end
time2 = toc(t1)

E = E/E(1);
w = [E,fliplr(conj(E))];
plot_evans(w)
winding_number(w)

E2 = E;

%--------------------------------------------------------------------------
% Evans function with continuous orthogonaliation
%--------------------------------------------------------------------------

E = zeros(size(preimage));


t1 = tic;
parfor j=1:length(E)
    E(j) = evans_co(s,p,preimage(j),e,m);
end
time3 = toc(t1)

E = E/E(1);
w = [E,fliplr(conj(E))];
plot_evans(w)
winding_number(w)

E3 = E;



figure;
hold on;
plot(E1,'.-k');
plot(E2,'o-r');
plot(E3,'*-g');


















