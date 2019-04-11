%% Set up variables for loading the profile.
addpath('../../Core');
clc; clear all; close all; beep off; curr_dir = cd;
orig_dir = cd;
folder='rns_code';
subfolder= 'data' ;
mathcalE = 7; % unstable wave
epsilon = 0.1;

%% Load the profile and create profile_fun, a deval for the profile.
file_name = ['profile_E_',num2str(10000000000*mathcalE),'eps',num2str(10000000000*epsilon)];
ld = load(['data/',file_name]);
d = ld.var;
p = d{1};
s = ld.var{2};
profile_fun = @(x) deval_loaded(x,s,p);

%% Test the profile
L = 10; % spatial infinity
num_homological_equations = 20;
scl = 0.2;

dom = linspace(-L,L,1600);

U_n = profile_fun(dom);

plot(dom, U_n);
return;

%% Parameterization parameters.
L = 10; % spatial infinity
num_homological_equations = 20;
scl = 0.2;
p.cnu = 1;
p.L = L;
p.Ti_weight = 0.99;

%% Test the profile is correct FIXME 

% x = linspace(-L,L,8001);
% 
% u = 1-3*p.gamma./(1+p.Q*cosh(x/sqrt(p.gamma)));
% v =3./(1+p.Q*cosh(x/sqrt(p.gamma)));
% 
% 
% delx = x(2)-x(1);
% u_xx = (u(3:end)-2*u(2:end-1)+u(1:end-2))/(delx^2);
% v_xx = (v(3:end)-2*v(2:end-1)+v(1:end-2))/(delx^2);
% 
% u = u(2:end-1);
% v = v(2:end-1);
% eq1 = u_xx-u.*v.^2+p.alpha*(1-u);
% eq2 = v_xx+u.*v.^2/p.gamma-v/p.gamma;
% 
% max(abs(eq1))
% max(abs(eq2))
% return;


%% Test profile_fun is correct FIXME

% x = linspace(-L,L, 100);
% eps = 1e-10;
% delX = x(2)-x(1);
% profile_fun = @(x)[1-3*p.gamma./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
%     -3*p.Q*sqrt(p.gamma)*sinh(x/sqrt(p.gamma))./(1+p.Q*cosh(x/sqrt(p.gamma))).^2;...
%     3./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
%     3*p.Q*(1/sqrt(p.gamma))*sinh(x/sqrt(p.gamma))./(1+p.Q*cosh(x/sqrt(p.gamma))).^2];
%     
% profile_fun(1)
% profileFunFd = (profile_fun(x+eps)-profile_fun(x-eps))./eps;
% uPrimeFd = profileFunFd(1);
% vPrimeFd = profileFunFd(3);
% profileFun = profile_fun(x);
% uPrime = profileFun(2);
% vPrime = profileFun(4);
% 
% max(abs(uPrimeFd - uPrime))
% max(abs(vPrimeFd - vPrime))
% 
% return;


%% solve for the eigenfunction
s.profile_fun = profile_fun;
s.profile_der = @(x, h)((profile_fun(x+h)-profile_fun(x-h))/h);
eigfun_ode = @(x,y,lambda)([A(x,lambda,s,p)*y(1:7); ...
    A(x+L,lambda,s,p)*y(8:14)]);
fun_bc = @(ya,yb,lambda)(bc_eigfun(ya,yb,lambda,p));

guess_fun = @(x)[sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];
lambda0 = 1;

guess = @(x)[guess_fun(x);guess_fun(x+L)];

solinit = bvpinit(linspace(-L,0,30),guess,lambda0);
options = bvpset('RelTol',1e-12,'AbsTol',1e-12,'Nmax', 20000);
eigfun_sol = bvp5c(eigfun_ode,fun_bc,solinit,options);
lam = eigfun_sol.parameters;


%% plot eigenfunction

hold on;
x = linspace(-L,0,500);
y = deval(eigfun_sol,x);
plot(x,y(1:4,:),'-k','LineWidth',2)
plot(x+L,y(5:8,:),'-k','LineWidth',2)
return;

%% check that eigenfunction is correct
% 
% x = linspace(-L,0,1000);
% delx = x(2)-x(1);
% u = zeros(length(x),1);
% v = u;
% u_der = u;
% v_der = u;
% uhat = u;
% vhat = u;
% u_der_der = u;
% v_der_der = u;
% for j = 1:length(x)
%    [temp,temp_der] = deval(eigfun_sol,x(j));
%    u(j) = temp(1);
%    v(j) = temp(3);
%    u_der(j) = temp(2);
%    v_der(j) = temp(4);
%    u_der_der(j) = temp_der(2);
%    v_der_der(j) = temp_der(4);
%    uhat(j) = 1-3*p.gamma./(1+p.Q*cosh(x(j)/sqrt(p.gamma)));
%    vhat(j) =3./(1+p.Q*cosh(x(j)/sqrt(p.gamma)));
% end
% 
% eq1 = -lam*u+u_der_der-u.*vhat.^2-2*uhat.*vhat.*v-p.alpha*u;
% eq2 = -lam*v+v_der_der+vhat.^2.*u/p.gamma+2*uhat.*vhat.*v/p.gamma-v/p.gamma;
% max(abs(eq1))
% max(abs(eq2))
% return 


%% Chebyshev interpolation

% homological equations to be solved
num_homological_equations = 1; %10;
%scl = 0.2; 
scl = 0.2;
n_vals = 2:1:num_homological_equations;

% degree of polynomial
N = 1001;
a = -L;
b = L;

% Create Chebyshev nodes
theta = ((0:1:N-1)+0.5)*pi/N;
x0 = cos(theta);
x = 0.5*(a+b)+0.5*(b-a)*x0; 

% Transformation to get Chebyshev coefficients
Id2 = (2/N)*speye(N);
Id2(1,1) = Id2(1,1)/2;
Tcf = Id2*cos(theta.'*(0:1:N-1)).';

profile_fun = @(x) extend(profileDeval, x, -L, L);

%% For P=1, get the chebyshev and p values.
fprintf("P = 1\n");
F0 = zeros(N,4);
for j = 1:N
   F0(j,:) = profile_fun(x(j));
end

%plot(x,F0,'-k','MarkerSize',18)

for j = 1:4
    [cf,fun] = get_chebyshev_coefficients(F0(:,j),a,b,'first kind');
    P{1}{j} = fun;
%    max(max(abs(F0(:,j)-fun(x))))
%    figure
%    hold on
%    plot(x,F0(:,j),'-k');
   % plot(linspace(a,b,1000),fun(linspace(a,b,1000)),'-r');
end
% return


%% For P=2, get the chebyshev and p values.
fprintf("P = 2\n");
F1 = zeros(N,4);
for j = 1:N 
    if x(j) <= 0
        temp = deval(eigfun_sol,-abs(x(j)));
        F1(j,:) = temp(1:4);
    else
        temp = deval(eigfun_sol,-L+abs(x(j)));
        F1(j,:) = temp(5:8);
    end
end

%F1(1:1:N,:) = F1(N:-1:1,:); % flip x.
%F1(1:1:N,:) = -F1(1:1:N,:); %flip y.

for j = 1:4
    [cf,fun] = get_chebyshev_coefficients(scl*F1(:,j),a,b,'first kind');
    P{2}{j} = fun;
%     max(max(abs(scl*F1(:,j)-fun(x))))
%     figure
%     hold on
%     plot(x,F1(:,j),'-k');
%     plot(x,fun(x),'-r');
end

% return


%% For P=3...N, get the chebyshev and p values.

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);
cnt = 0;
for n = n_vals
    
    
    fprintf('P =%4g\n',n+1);
    
    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam,p));
    bc_handle = @(ya,yb)(bc(ya,yb,n,lam,p)); 
    total_fun = @(x)[P{n}{1}(x);P{n}{2}(x);P{n}{3}(x);P{n}{4}(x)];
    
    solinit = bvpinit(linspace(a,b,30),total_fun);
    sol = bvp5c(ode_handle,bc_handle,solinit,options);
        
    temp = deval(sol,x);
   
    for j = 1:4
        [cf,fun] = get_chebyshev_coefficients(temp(j,:),a,b,'first kind');
        P{n+1}{j} = fun;
%         figure
%         hold on
%         plot(x,temp,'-k');
%         drawnow;
    end
    
%     return
end


%% Create y0

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

%% Test the values of P = 2

% sigval = 1;
% xVals = linspace(a,b,20*round(b-a)+1);
% y1 = P{2}{1}(xVals).';
% y2 = P{2}{3}(xVals).';
% u  = zeros(2,length(xVals));
% u = u+sigval*[y1;y2];
% plot(xVals, u);
% return;

%% Get u0
%
% test the method
%
%figure;
%hold on;

sigma0 = 0.1 %0.1
xgrid = linspace(a,b,10*round(b-a)+1);
u0 = zeros(2,length(xgrid));
for j = 0:num_homological_equations
    y1 = P{j+1}{1}(xgrid).';
    y2 = P{j+1}{3}(xgrid).';
    u0 = u0+sigma0^j*[y1;y2]; %
end
% plot(xgrid,u0,'-r','LineWidth',2);
% return

% plot(xgrid,u0,'-k','LineWidth',2);
% plot(-fliplr(xgrid),fliplr(u0),'-k','LineWidth',2);


%% Get u1 through interpolation.
sigma1 = 0.2 %0.2
u1 = zeros(2,length(xgrid));
for j = 0:num_homological_equations
    y1 = P{j+1}{1}(xgrid).';
    y2 = P{j+1}{3}(xgrid).';
    u1 = u1+sigma1^j*[y1;y2]; 
end
%plot(xgrid,u1,'--r','LineWidth',2); %
%plot(-fliplr(xgrid),fliplr(u1),'--r','LineWidth',2); %

% figure;
% hold on;
% plot(xgrid,u0,'-r','LineWidth',2);
% plot(xgrid,u1,'-b','LineWidth',2);
% plot(xgrid,P{1}{1}(xgrid),'--m');
% axis([-5,5,-0.1,1.6]);
% drawnow;
% return


%% TIME EVOLUTION
% xgrid = [xgrid,-fliplr(xgrid(1:end-1))];
bc_L = [P{1}{1}(-L);P{1}{3}(-L)];
bc_R = bc_L;

bc_L_fun = @(U_n,U_o,K,H,p)(U_n(:,1)-bc_L);
bc_R_fun = @(U_n,U_o,K,H,p)(U_n(:,end)-bc_R);
bc_L_jac_fun = @(U_n,U_o,K,H,p)eye(2);
bc_R_jac_fun = @(U_n,U_o,K,H,p)eye(2);

T = (log(sigma1)-log(sigma0))/lam

TPoints = 30;
% return

tgrid = linspace(0,T,TPoints);
K = tgrid(2)-tgrid(1);
H = xgrid(2)-xgrid(1);

U0 = u0;
U1 = u1;

% U0 = [u0,fliplr(u0(:,1:end-1))];
% U1 = [u1,fliplr(u1(:,1:end-1))];

U_n = U0;
U_o = U0;

tol = 1e-18;
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
    % plot(xgrid,prof,'-g','LineWidth',2);
    % axis([-.25,.25,1.65,1.9]);    
    % axis([-5,5,-0.1,1.6]);
    axis([-1,1,-0.1,1.9]);
    drawnow;
    pause(0.1);

    U_o = U_n;

end






