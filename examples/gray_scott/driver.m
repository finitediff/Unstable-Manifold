clc; clear all; close all; beep off; curr_dir = cd;
addpath(strcat(pwd,'/../../core/matlab'))

%% User defined variables

% System variables
p.gamma = 1/9; % 0 < gamma < 2/9 

% Controls
L = 10;
scl = 0.2;
num_homological_equations = 20;

%dependent variables (alpha*gamma = 1)
p.alpha = 1/p.gamma;
p.Q = sqrt(1-9*p.gamma/2);

%% solve for the eigenfunction
fprintf("n=0, Loading profile\n\n");
profile_fun = @(x)[1-3*p.gamma./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
    -3*p.Q*sqrt(p.gamma)*sinh(x/sqrt(p.gamma))/(1+p.Q*cosh(x/sqrt(p.gamma)))^2;...
    3./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
    3*p.Q*(1/sqrt(p.gamma))*sinh(x/sqrt(p.gamma))/(1+p.Q*cosh(x/sqrt(p.gamma)))^2];
    
fprintf("n=1, Loading eigenfunction\n\n");
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

%% plot eigenfunction
% x = linspace(-L,0,500);
% y = deval(eigfun_sol,x);
% hold on;
% plot(x,y(1:4,:),'-g','LineWidth',2)
% plot(x+L,y(5:8,:),'-b','LineWidth',2)

% return


%% check that eigenfunction is correct
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
% 
% max(abs(eq1))
% max(abs(eq2))
% 
% return


%% system info

% n = 0;
% 
% B = [0, 1, 0, 0; 
%     (lam*n+p.alpha), 0, 0, 0;
%     0, 0, 0, 1;
%     0, 0, lam*n+1/p.gamma, 0];
% 
% [V,D] = eig(B)


%% program guts

% homological equations to be solved
n_vals = 2:1:num_homological_equations;

% degree of polynomial
N = 1001;
a = -L;
b = L;

% Chebyshev nodes in [a,b]
theta = ((0:1:N-1)+0.5)*pi/N;
x0 = cos(theta);
x = 0.5*(a+b)+0.5*(b-a)*x0; 

% Transformation to get Chebyshev coefficients
Id2 = (2/N)*speye(N);
Id2(1,1) = Id2(1,1)/2;
Tcf = Id2*cos(theta.'*(0:1:N-1)).';

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

% Get the Chebyshev coefficients for the eigenfunction (on the interal [-L,0])
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

% x = linspace(-L,0,500);
% y = deval(eigfun_sol,x);
% hold on;
% plot(x,scl*y(1:4,:),'--g','LineWidth',2)
% plot(x+L,scl*y(5:8,:),'--b','LineWidth',2)
% return

% figure;
% hold on;
% sigma0 = 0.1 % theta*exp(lambda*t)
% u0 = zeros(2,length(x));
% for j = 1%:num_homological_equations
%     y1 = P{j+1}{1}(x).';
%     y2 = P{j+1}{3}(x).';
%     u0 = u0+sigma0^j*[y1;y2];
% end
% % plot(x,u0,'-r','LineWidth',2);


options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);

cnt = 0;
for n = n_vals
    fprintf(strcat('n=', num2str(n), ', Generating function\n\n'));
    
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


%% Finite difference test
fprintf('Preparing graphs\n\n');
figure;
hold on;

%% User defined FD parameters
t1 = 0.25;
xpoints = 801; % 201 401 801
tpoints = 1601; % 101 401 1601
xgrid = linspace(-10,10, xpoints);
tgrid = linspace(0,t1,tpoints);
sigma0 = 0.1;
sigma1 = sigma0*exp(lam*t1);

%% Generate u0 and u1.

u0 = zeros(2,length(xgrid));
u1 = zeros(2,length(xgrid));
for j = 0:num_homological_equations
    y1 = P{j+1}{1}(xgrid).';
    y2 = P{j+1}{3}(xgrid).';
    u0 = u0+sigma0^j*[y1;y2];
    u1 = u1+sigma1^j*[y1;y2];
end

% Create finite difference boundary conditions
bc_L = [P{1}{1}(-L);P{1}{3}(-L)];
bc_R = [P{1}{1}(L);P{1}{3}(L)];
bc_L_fun = @(U_n,U_o,K,H,p)(U_n(:,1)-bc_L);
bc_R_fun = @(U_n,U_o,K,H,p)(U_n(:,end)-bc_R);
bc_L_jac_fun = @(U_n,U_o,K,H,p)eye(2);
bc_R_jac_fun = @(U_n,U_o,K,H,p)eye(2);


T = (log(sigma1)-log(sigma0))/lam;
K = tgrid(2)-tgrid(1);
H = xgrid(2)-xgrid(1);

U0 = u0;
U1 = u1;
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
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    clf;
    hold on;
    plot(xgrid,U0,'-r','LineWidth',2);
    plot(xgrid,U1,'-b','LineWidth',2);
    plot(xgrid,U_n,'--k','LineWidth',2);
    axis([-L,L,-0.1,2]);
    drawnow;
    pause(0.1);

    U_o = U_n;

end






