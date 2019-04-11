clc; clear all; close all; beep off; curr_dir = cd;

% Nagumo - reaction diffusion equation
% 
% u_t = u_{xx} - u + u^3

% 
% controls
%

L = 20;
scl = 0.1;
num_homological_equations = 30;

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

P{1}.fun = @(x)([sqrt(2)*sech(x); ...
                sqrt(2)*(-tanh(x).*sech(x)); ...
                sqrt(2)*sech(-x); ...
                sqrt(2)*(-tanh(-x).*sech(-x))
                ]);

% eigenfunction
P{2}.fun = @(x)([scl*eig_fun(x);scl*eig_fun(-x)]);

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);

cnt = 0;
for n = n_vals
    
    fprintf('n = %4g\n',n);
    
    ode_handle = @(x,y)(p_ode_double(x,y,P,n,lam));

    bc_handle = @(ya,yb)(bc_double(ya,yb,n)); 

    solinit = bvpinit([-L,0],P{n}.fun);

    sol = bvp5c(ode_handle,bc_handle,solinit);
        
    P{n+1}.mu_L = sqrt(3*n+1);
    P{n+1}.mu_R = -sqrt(3*n+1);
    P{n+1}.sol = sol;
    
    % degree of polynomial
    N = 120;
    
    a = -L;
    b = 0;

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
    cf = Tcf*interp_double(sol,x,n).';
    
    P{n+1}.fun = @(x)(interp_cheby(cf,x,n,N,a,b));

end

x = linspace(-L,0);

hold on;
for j = 1:num_homological_equations
    fprintf('j = %4g\n',j);
    y = P{j}.fun(x);
    plot(x,y(1:2,:),-x,y(3:4,:))
    nrm = max(max(abs(y)))
end

%
% test the method
%

figure;
hold on;

%sigma = sigma ;
xgrid = linspace(-L,0,10*L+1);
u0 = zeros(1,2*length(xgrid)-1);
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u0 = u0+sigma^j*[y(1,:),fliplr(y(3,2:end))];
end
xgrid = [xgrid,-fliplr(xgrid(2:end))];
plot(xgrid,u0,'-k');


sigma = 1;
u1 = zeros(1,2*length(xgrid)-1);
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u1 = u1+sigma^j*[y(1,:),fliplr(y(3,2:end))];
end
xgrid = [xgrid,-fliplr(xgrid(2:end))];
plot(xgrid,u1,'-g');


T = log(2)/3;


tgrid = linspace(0,T,131);
delt = tgrid(2)-tgrid(1);
delx = xgrid(2)-xgrid(1);

M = zeros(length(tgrid),length(xgrid));

M(1,:) = u0;
for j = 1:length(tgrid)-1
   for k = j+1:length(xgrid)-j
      M(j+1,k) = M(j,k) +(delt/delx^2)*(M(j,k+1)-2*M(j,k)+M(j,k-1))+ ...
          delt*(M(j,k)+M(j,k)^3);
   end
end

plot(xgrid,M(end,:),'-b');



