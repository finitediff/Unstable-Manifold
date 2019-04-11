
%% Nagumo - reaction diffusion equation ( U_t = U_{xx} - U + U^3 )

% A driver that finds the unstable manifold of Nagumo's equation found by
% Leaving negative infinity in the direction of the eigenfunction.
% Compares the results with the finite difference method found in stablab
% in order to confirm accuracy.  To run this file, the user must:
%
% 1. Add finite_diff_advance, a function found in the time evolution 
% portion of stablab.
%
% 2. Add the "core/matlab" folder found in this package.  This folder
% contains functions needed to solve chebychev polynomials.

clc; clear all; close all; beep off; curr_dir = cd;

%% User defined parameters

L = 20; % Spatial infinity.  
scl = 0.3; % Scale of eigenfunction.
num_homological_equations = 30; % Number of homological equations.

%% Set up variables

% System specific variables
lam = 3; % eigenvalue of linearized PDE

% Other variables
n_vals = 2:1:num_homological_equations;  % homological equations.

%% Find the homological equations, P{n}.fun

% Unstable wave solution with coordinates to interval [-L,0]
P{1}.fun = @(x)([sqrt(2)*sech(x); ...
                sqrt(2)*(-tanh(x).*sech(x))]);

% Eigenfunction
P{2}.fun = @(x) scl*eig_fun(x);

% Solve for each P{n} recursively
options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);
cnt = 0;
for n = n_vals
    
    fprintf('n = %4g\n',n);
    
    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam));
    bc_handle = @(ya,yb)(bc_fun(ya,yb,n)); 

    solinit = bvpinit([-L,L],P{n}.fun);
    sol = bvp5c(ode_handle,bc_handle,solinit);
        
    P{n+1}.mu_L = sqrt(3*n+1);
    P{n+1}.mu_R = -sqrt(3*n+1);
    P{n+1}.sol = sol;
    
    % degree of polynomial
    N = 120;
    
    a = -L;
    b = L;

    % Chebyshev nodes
    theta = ((0:1:N-1)+0.5)*pi/N;
    x0 = cos(theta); % in [-1,1]
    x = 0.5*(a+b)+0.5*(a-b)*x0; % in [a,b]

    % Transformation to get Chebyshev coefficients
    Id2 = (2/N)*speye(N);
    Id2(1,1) = Id2(1,1)/2;
    Tcf = Id2*cos(theta.'*(0:1:N-1)).';
    cf = Tcf*deval(sol,x).';
    
    P{n+1}.fun = @(x)(interp_cheby_fun(cf,x,n,N,a,b));

end

%% Plot the homological equations

x = linspace(-L,L,1000);

hold on;
for j = 1:num_homological_equations
    fprintf('j = %4g\n',j);
    y = P{j}.fun(x);
    plot(x,y,'LineWidth',2)
    nrm = max(max(abs(y)))
end


%% Test the method

figure;
hold on;

% Plot the solution at sigma0
sigma0 = 0.1;
xgrid = linspace(-L,L,40*L+1);
u0 = zeros(1,length(xgrid));
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u0 = u0+sigma0^j*y(1,:);
end
plot(xgrid,u0,'-k','LineWidth',2);

% Plot the solution at sigma1
sigma1 = 0.7;
u1 = zeros(1,length(xgrid));
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u1 = u1+sigma1^j*y(1,:);
end
plot(xgrid,u1,'-g','LineWidth',2);

%% Time evolution

% Set up boundary conditions.
bc_L = sqrt(2) * sech(xgrid(1));
bc_R = sqrt(2) * sech(xgrid(end));
bc_L_fun = @(U_n,U_o,K,H,p) (U_n(1,1)-bc_L);
bc_R_fun = @(U_n,U_o,K,H,p) (U_n(1,end)-bc_R);
bc_L_jac_fun = @(U_n,U_o,K,H,p) 1;
bc_R_jac_fun = @(U_n,U_o,K,H,p) 1;

% Set up grid
T = (log(sigma1)-log(sigma0))/3
tgrid = linspace(0,T,200);
K = tgrid(2)-tgrid(1);
H = xgrid(2)-xgrid(1);


U_n = u0;
U_o = u0;
tol = 1e-8;
p = struct;

for j = 1:length(tgrid)
    
    % Update finite difference
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    % Plot current state
    j
    clf;
    hold on;
    plot(xgrid,u1,'-g','LineWidth',2);
    plot(xgrid,u0,'-k','LineWidth',2);
    plot(xgrid,U_n,'--r','LineWidth',2);
    drawnow;
    pause(0.1);
    U_o = U_n;

end













