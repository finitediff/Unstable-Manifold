clc; clear all; close all; beep off; curr_dir = dir;

%% User defined parameters.

% Parameters
p.a = 3; % a > 0
p.gamma = 7; % gamma >= 1

% Controls
L = 299; %300
num_homological_equations = 12;
bvp_points = 30;
N = 0; % degree of cheb polynomial. Zero for deval method.

%% Get all parameters
fprintf(strcat('Manifold for Shrodingers equation with parameters, \n', ...
    'a=', num2str(p.a,'%.1f'), ', gamma=', num2str(p.gamma,'%.1f'), ', L=', ...
    num2str(L, '%.1f'), ', ChebNodes(N)=', num2str(N), ', HomEqs=', ...
    num2str(num_homological_equations),'\n\n'));

scl = 0.2;
p.mu = (p.a-sqrt(p.gamma^2-1))/(p.a+sqrt(p.gamma^2-1));
p.nu = 1/(p.a+sqrt(p.gamma^2-1));

%% Set up profile solution as a function.
fprintf('n=0, Loading function\n\n');
profile_fun = @(x) profile(x);

%% Load the eigenvalue and eigenfunction.
fprintf('n=1, Loading function\n\n');
[eig_fun_left, lam] = eigenfunction();
eig_fun = @(x) deval_left(eig_fun_left, x).';

% homological equations to be solved
n_vals = 2:1:num_homological_equations;

%% Chebyshev interpolation
a = -L;
b = L;

%% Get the Chebyshev coefficients for the profile (on the interal [-L,L])
P{1} = profile_fun;
P{2} = eig_fun;

%% Find the homological equations, loading, generating, saving

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);
cnt = 0;

for n = n_vals
    
    % data_a_gamma_L_chebychev_homologicalEq.mat
    fileName = strcat(strrep(strcat('data-', num2str(p.a,'%.10f'), '-', ...
        num2str(p.gamma,'%.10f'), '-', num2str(L, '%.2f'), '-', ...
        num2str(N), '-',num2str(bvp_points),'-',num2str(n)), '.', 'p'), '.mat');
    
    % If it exists,
    if exist(fileName, 'file') == 2
        fprintf(strcat('n=', num2str(n), ', Loading function\n\n'));
        ld = load(fileName);
        P{n+1} = ld.homEq;
        
    % If it doesn't exist
    else
        fprintf(strcat('n=', num2str(n), ', Generating function\n\n'));
        ode_handle = @(x,y)(ode_fun_deval(x,y,P,n,lam,p));
        bc_handle = @(ya,yb)(bc(ya,yb,n,lam,p)); 
        total_fun = @(x)[P{n}(x)];
        solinit = bvpinit(linspace(-L,L,bvp_points),total_fun);
        sol = bvp5c(ode_handle,bc_handle,solinit,options);  

        P{n+1} = @(x) deval(sol, x);
        homEq = P{n+1};
        save(fileName, 'homEq');
    end
    

end



%% User defined parameters for finite difference
xgrid = linspace(-1,1, 1000);
tpoints = 200;
sigma0 = 0.6;
sigma1 = 0.61;

%% Dependant parameters

% Prepare u0 and u1
fprintf('Preparing graphs\n\n');
figure;
hold on;
u0 = zeros(2, length(xgrid));
u1 = zeros(2, length(xgrid));

for j = 0:num_homological_equations;
    temp = P{j+1}(xgrid);
    y1 = temp(1,:);
    y2 = temp(3,:);
    u0 = u0 + sigma0^j*[y1;y2];
    u1 = u1 + sigma1^j*[y1;y2];
end
   
% Plot u0 and u1
plot(xgrid, real(u0));
plot(xgrid, real(u1));
drawnow;

%% Prepare for finite difference.
T = (log(sigma1)-log(sigma0))/real(lam)
tgrid = linspace(0,T,tpoints);
lambdagrid = exp((tgrid)*real(lam)+log(sigma0));
K = tgrid(1,2) - tgrid(1,1);
H = xgrid(1,2) - xgrid(1,1);
tol = 1e-10;
un = u0;


%% Plot finite difference progression.
for i = 2:tpoints
    i
    
    % Update the boundary values.
    sigman = lambdagrid(i);
    bc_L = [0;0];
    bc_R = [0;0];
    d_bc_L = [0;0];
    d_bc_R = [0;0];
    
    % Produce the boundary conditions.
    for j = 0:num_homological_equations;
        temp_L = P{j+1}(xgrid(1,1));
        temp_R = P{j+1}(xgrid(1,end));
    
        bc_L = bc_L + sigman^j*[temp_L(1,:); temp_L(3,:)];
        bc_R = bc_R + sigman^j*[temp_R(1,:); temp_R(3,:)];
        d_bc_L = d_bc_L + sigman^j*[temp_L(2,:); temp_L(4,:)];
        d_bc_R = d_bc_R + sigman^j*[temp_R(2,:); temp_R(4,:)];
    end
    
    % Create boundary functions
    bc_L_fun = @(U_n,U_o,K,H,p)(U_n(:,1)-bc_L);
    bc_R_fun = @(U_n,U_o,K,H,p)(U_n(:,end)-bc_R);
    bc_L_jac_fun = @(U_n,U_o,K,H,p)eye(2);
    bc_R_jac_fun = @(U_n,U_o,K,H,p)eye(2);
    
    % Update the initial guess.
    uo = un;
    un = finite_diff_advance(un,uo,K,H,p,tol,@fd_F,@fd_jac,bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);
    
    clf;
    hold on;
    plot(xgrid, real(u0));
    plot(xgrid, real(u1));
    plot(xgrid, real(un));
    drawnow;
    %pause(0.01);
end

% I have checked the boundary conditions, and I've checked the inputted PDE
% I found a small error in the FD code for expressing the U term in crank-
% nicholson.  This has been fixed.  
% Next I should check
