clc; clear all; close all; beep off; curr_dir = cd;

% Nagumo - reaction diffusion equation
% 
% u_t = u_{xx} - u + u^3

% 
% controls
%

L = 20;
scl = 0.3;
num_homological_equations = 40;

% 
% program guts
%

% eigenvalue of the linearized PDE
lam = 3;

% homological equations to be solved
n_vals = 2:1:num_homological_equations;

% left and right spatial infinity
XL = -L;
XR = L;

% change of coordinates to interval [-L,0]
P{1}.fun = @(x)([sqrt(2)*sech(x); ...
                sqrt(2)*(-tanh(x).*sech(x))]);

% eigenfunction
P{2}.fun = @(x)scl*eig_fun(x);

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
    
    % nodes in [-1,1] 
    x0 = cos(theta);

    % nodes in [a,b]
    x = 0.5*(a+b)+0.5*(a-b)*x0; 

    % Transformation to get Chebyshev coefficients
    Id2 = (2/N)*speye(N);
    Id2(1,1) = Id2(1,1)/2;
    Tcf = Id2*cos(theta.'*(0:1:N-1)).';
    cf = Tcf*deval(sol,x).';
    
    P{n+1}.fun = @(x)(interp_cheby_fun(cf,x,n,N,a,b));

end

x = linspace(-L,L);

hold on;
for j = 1:num_homological_equations
    fprintf('j = %4g\n',j);
    y = P{j}.fun(x);
    plot(x,y)
    nrm = max(max(abs(y)))
end


% return
%
% test the method
%

figure;
hold on;

sigma0 = 0.1; % sigma e ^ (lambda theta)
xgrid = linspace(-L,L,40*L+1);
u0 = zeros(1,length(xgrid));
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u0 = u0+sigma0^j*y(1,:);
end
plot(xgrid,u0,'-k','LineWidth',2);


sigma1 = 0.7;
u1 = zeros(1,length(xgrid));
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u1 = u1+sigma1^j*y(1,:);
end
plot(xgrid,u1,'-g','LineWidth',2);

drawnow;

% theta = sigma e^(lambda theta)
% theta ^ n
% solution = p_n theta^n

%Create the solution.

for curr_num_homological_equations = 0:1:num_homological_equations
xgrid = linspace(-2,2,40*L+1);
thetaMax = 0.9;
numPoints = 100;
solution = zeros(1,length(xgrid));
for i=0:1:numPoints
    theta = thetaMax*i/numPoints;
    
    u0 = zeros(1,length(xgrid));
    
    for j = 0:curr_num_homological_equations
        y = P{j+1}.fun(xgrid);
        u0 = u0+theta^j*y(1,:);
    end
    
    solution(i+1,:) = u0;
end

    dlmwrite(strcat('temp_',int2str(curr_num_homological_equations),'.txt'), solution)
end
