clc; clear all; close all; beep off; curr_dir = dir;

%% User defined parameters.

% Parameters
p.a = 3; % a > 0
p.gamma = 7; % gamma >= 1

% Controls
L = 300;
num_homological_equations = 12;
bvp_points=2000;
N = 0; % degree of cheb polynomial

%% Get all parameters
fprintf(strcat("Manifold for Shrodinger's equation with parameters, \n", ...
    "a=", num2str(p.a,'%.1f'), ', gamma=', num2str(p.gamma,'%.1f'), ", L=", ...
    num2str(L, '%.1f'), ', ChebNodes(N)=', num2str(N), ", HomEqs=", ...
    num2str(num_homological_equations),"\n\n"));

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

%% test the method

% Plot
fprintf("Preparing graphs\n\n");
figure;
hold on;

% Choose sigma values
sigma = linspace(0,0.5,2);
xgrid = linspace(-300,300,10000);

for i = 1:size(sigma,2);
    i;
    currSigma = sigma(1,i);
    u = zeros(2,length(xgrid));
    for j = 0:num_homological_equations;
        temp = P{j+1}(xgrid);
        y1 = temp(1,:);
        y2 = temp(3,:);
        u = u+currSigma^j*[y1;y2];
    end
    
    plot(xgrid, real(u));
end

drawnow;

% Format for saving files is,
% data_a_gamma_L_chebychev_homologicalEq


