clc; clear all; close all; beep off; curr_dir = cd;
%
%%% This driver was directly copied from the Gray scott example and
%%% modified to have the functions from the Shrodinger system.  
%


p.a = 3; % a > 0
p.gamma = 7; % gamma >= 1
p.Q = 1; % Delete This.
p.alpha = p.a;

%
% dependent variables
%

p.mu = (p.a-sqrt(p.gamma^2-1))/(p.a+sqrt(p.gamma^2-1));
p.nu = 1/(p.a+sqrt(p.gamma^2-1));

% 
% controls
%

L = 10; % spatial infinity
scl = 0.2;
num_homological_equations = 20;

% setting a temporary eigenfunction.  Replace this with
% The correct eigenfunction.
s = struct;
eigfun_ode = @(x,y,lambda)([A_explicit(x,lambda,s,p)*y(1:4); ...
    A_explicit(x+L,lambda,s,p)*y(5:8)]);
fun_bc = @(ya,yb,lambda)(bc_eigfun(ya,yb,lambda,p));

guess_fun = @(x)[sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];
lambda0 = 1;

guess = @(x)[guess_fun(x);guess_fun(x+L)];

solinit = bvpinit(linspace(-L,0,30),guess,lambda0);
options = bvpset('RelTol',1e-12,'AbsTol',1e-12,'Nmax', 20000);
eigfun_sol = bvp5c(eigfun_ode,fun_bc,solinit,options);
lam = eigfun_sol.parameters;

%% program guts
% 
% 
%

% homological equations to be solved
n_vals = 2:1:num_homological_equations;


%% Set up Chebyshev interpolation

% degree of polynomial
N = 1001;

a = -L;
b = L;

% Chebyshev nodes
theta = ((0:1:N-1)+0.5)*pi/N;

% nodes in [-1,1] 
x0 = cos(theta);

% nodes in [a,b]
x = 0.5*(a+b)+0.5*(b-a)*x0; 

% Transformation to get Chebyshev coefficients
Id2 = (2/N)*speye(N);
Id2(1,1) = Id2(1,1)/2;
Tcf = Id2*cos(theta.'*(0:1:N-1)).';

%% Get chebychev coefficients for profile

% Set a temporary profile solution.
profile_fun = @(x)[1-3*p.gamma./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
    -3*p.Q*sqrt(p.gamma)*sinh(x/sqrt(p.gamma))/(1+p.Q*cosh(x/sqrt(p.gamma)))^2;...
    3./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
    3*p.Q*(1/sqrt(p.gamma))*sinh(x/sqrt(p.gamma))/(1+p.Q*cosh(x/sqrt(p.gamma)))^2];

% Get the Chebyshev coefficients for the profile (on the interal [-L,L])
F0 = zeros(N,4);
for j = 1:N
   F0(j,:) = profile_fun(x(j));
end
% plot(x,F0,'-k','MarkerSize',18)
% return

for j = 1:4
    [cf,fun] = get_chebyshev_coefficients(F0(:,j),a,b,'first kind');
    P{1}{j} = fun;
%     max(max(abs(F0(:,j)-fun(x))))
%     figure
%     hold on
%     plot(x,F0(:,j),'-k');
%     plot(linspace(a,b,1000),fun(linspace(a,b,1000)),'-r');
end

% figure;
% x = linspace(-L,L,1001);
% plot(x,P{1}{3}(x),'-k');

%% Get the Chebyshev coefficients for the eigenfunction (on the interal [-L,0])
F1 = zeros(N,4);
for j = 1:N
    
    if x(j) <= 0
        temp = deval(eigfun_sol,x(j));
        F1(j,:) = temp(1:4);
    else
        temp = deval(eigfun_sol,x(j)-L);
        F1(j,:) = temp(5:8);
    end
end

for j = 1:4
    [cf,fun] = get_chebyshev_coefficients(scl*F1(:,j),a,b,'first kind');
    P{2}{j} = fun;
%     max(max(abs(scl*F1(:,j)-fun(x))))
%     figure
%     hold on
%     plot(x,F1(:,j),'-k');
%     plot(x,fun(x),'-r');
end

%% Get chebychev coefficients for homological equations.

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);

cnt = 0;
for n = n_vals
    
    fprintf('n = %4g\n',n);
    
    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam,p));
    
    bc_handle = @(ya,yb)(bc(ya,yb,n,lam,p)); 

    total_fun = @(x)[P{n}{1}(x);P{n}{2}(x);P{n}{3}(x);P{n}{4}(x)];
    
    solinit = bvpinit(linspace(a,b,30),total_fun);

    sol = bvp5c(ode_handle,bc_handle,solinit,options);
        
    temp = deval(sol,x);
   
    
    for j = 1:4
        [cf,fun] = get_chebyshev_coefficients(temp(j,:),a,b,'first kind');
        P{n+1}{j} = fun;

    end
end



% x = linspace(a,b);

hold on;
y = zeros(4,length(x));
for j = 1:num_homological_equations
    fprintf('j = %4g\n',j);
    for k = 1:4
        y(k,:) = P{j}{k}(x);
    end
%     size(y)
%     plot(x,y)
    nrm = max(max(abs(y)))
end



%
% test the method
%

figure;
hold on;

sigma0 = 0.1 % theta*exp(lambda*t)
xgrid = linspace(a,b,20*round(b-a)+1);
u0 = zeros(2,length(xgrid));
for j = 0:num_homological_equations
    y1 = P{j+1}{1}(xgrid).';
    y2 = P{j+1}{3}(xgrid).';
    u0 = u0+sigma0^j*[y1;y2];
end
% plot(xgrid,u0,'-r','LineWidth',2);
% return

% plot(xgrid,u0,'-k','LineWidth',2);
% plot(-fliplr(xgrid),fliplr(u0),'-k','LineWidth',2);


sigma1 = 0.5
u1 = zeros(2,length(xgrid));
for j = 0:num_homological_equations
    y1 = P{j+1}{1}(xgrid).';
    y2 = P{j+1}{3}(xgrid).';
    u1 = u1+sigma1^j*[y1;y2];
end
% plot(xgrid,u1,'--r','LineWidth',2);
% plot(-fliplr(xgrid),fliplr(u1),'--r','LineWidth',2);

% figure;
% hold on;
% plot(xgrid,u0,'-r','LineWidth',2);
% plot(xgrid,u1,'-b','LineWidth',2);
% plot(xgrid,P{1}{1}(xgrid),'--m');
% axis([-5,5,-0.1,1.6]);
% drawnow;

% return
% TIME EVOLUTION

% xgrid = [xgrid,-fliplr(xgrid(1:end-1))];

bc_L = [P{1}{1}(-L);P{1}{3}(-L)];
bc_R = [P{1}{1}(L);P{1}{3}(L)];

bc_L_fun = @(U_n,U_o,K,H,p)(U_n(:,1)-bc_L);
bc_R_fun = @(U_n,U_o,K,H,p)(U_n(:,end)-bc_R);
bc_L_jac_fun = @(U_n,U_o,K,H,p)eye(2);
bc_R_jac_fun = @(U_n,U_o,K,H,p)eye(2);


T = (log(sigma1)-log(sigma0))/lam

% return

tgrid = linspace(0,T,200);
K = tgrid(2)-tgrid(1);
H = xgrid(2)-xgrid(1);

U0 = u0;
U1 = u1;

% U0 = [u0,fliplr(u0(:,1:end-1))];
% U1 = [u1,fliplr(u1(:,1:end-1))];

U_n = U0;
U_o = U_n;


clf;
hold on;
plot(xgrid,U0,'-r','LineWidth',2);
plot(xgrid,U1,'-b','LineWidth',2);
plot(xgrid,U_n,'--k','LineWidth',2);
axis([-L,L,-0.1,2]);
drawnow;
pause(0.1);
    

tol = 1e-8;
p.something = 0;
for j = 1:length(tgrid)
    
    j
    
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    clf;
    hold on;
    plot(xgrid,U0,'-r','LineWidth',2);
    plot(xgrid,U1,'-b','LineWidth',2);
    plot(xgrid,U_n,'--k','LineWidth',2);
%     plot(xgrid,prof,'-g','LineWidth',2);
    axis([-L,L,-0.1,2]);
    drawnow;
    pause(0.1);

    U_o = U_n;

end






