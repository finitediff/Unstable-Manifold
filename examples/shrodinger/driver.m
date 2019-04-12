clc; clear all; close all; beep off; curr_dir = dir;
addpath(strcat(pwd,'/../../core/matlab'))

%% User defined parameters.

% System Parameters (a > 0, gamma >= 1)
p.a = 3;
p.gamma = 7;
p.mu = (p.a-sqrt(p.gamma^2-1))/(p.a+sqrt(p.gamma^2-1));
p.nu = 1/(p.a+sqrt(p.gamma^2-1));

% Controls
L = 20;
num_homological_equations = 7;
chebPoints = 4001; 
eigScale = 0.2;
bvpPoints = 30;

% Chebychev setup
theta = ((0:1:chebPoints-1)+0.5)*pi/chebPoints;
x0 = cos(theta);
a = -L;
b = L;
x = 0.5*(a+b)+0.5*(b-a)*x0; 

%% Form P

% Get the profile
fprintf('n=0, Loading function\n\n');
P{1} = @(x) profile(x);

% Get the eigenfunction
fprintf('n=1, Loading function\n\n');
[eig_fun_left, lam] = eigenfunction();
P{2} = @(x) deval_left(eig_fun_left, x).';

% Form the homological equations
nvals = 2:num_homological_equations;
for n=nvals
    fprintf(strcat('n=', num2str(n), ','));
    getHomEqFast = optimize(@getHomEq);
    P = getHomEqFast(x,P,n,L,@ode_fun,@bc,lam,p,bvpPoints);
    fprintf("\n\n");
end

%% User defined parameters for finite difference
xgrid = linspace(-10,10, 1000);
tpoints = 100;
sigma0 = 0.1;
sigma1 = 0.21;

%% Dependant parameters for finite difference

% Should I be dividing by abs(lam) or real(lam) or just lam?
% The P{3} term explodes at the end points.  This is either a problem with
% the projective boundary conditions or with the chebychev interpolation.

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
T = (log(sigma1)-log(sigma0))/abs(lam)
tgrid = linspace(0,T,tpoints);
lambdagrid = exp((tgrid)*abs(lam)+log(sigma0));
K = tgrid(1,2) - tgrid(1,1);
H = xgrid(1,2) - xgrid(1,1);
tol = 1e-10;
un = u0;


%% Plot finite difference progression.
for i = 2:tpoints
    % tgrid(i)
    
    % Update the boundary values.
    sigman = lambdagrid(i);
    bc_L = [0;0];
    bc_R = [0;0];
    d_bc_L = [0;0];
    d_bc_R = [0;0];
    
    % Produce the boundary conditions.
    temp_L = zeros(4,1);
    temp_R = zeros(4,1);
    for j = 0:num_homological_equations;
        temp_L = P{j+1}(xgrid(1));
        temp_R = P{j+1}(xgrid(end));
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
    
    % Verify the accuracy of the method. 
    % It looks like it is performing correctly with error very small.
    %     u = (un+uo)/2;
    %     ux = (u(:,3:end) - u(:,1:end-2))/(2*H);
    %     ut = (un-uo)/K;
    %     uxx = (u(:,3:end) - 2*u(:,2:end-1) + u(:,1:end-2))/(H^2);
    %     u = u(:,2:end-1);
    %     ut = ut(:,2:end-1);
    %     temp = ut(1,:) + uxx(2,:) - p.mu*u(2,:)+(u(1,:).^2 + u(2,:).^2).*u(2,:);
    %     temp2 = ut(2,:) - uxx(1,:) + u(1,:) - (u(1,:).^2 + u(2,:).^2).*u(1,:)+2*p.nu*u(2,:);
    % 
    %     fprintf("error");
    %     max(temp(1,:))
    %     max(temp2(1,:))
   
    clf;
    hold on;
    plot(xgrid, real(u0), 'LineWidth',2);
    plot(xgrid, real(u1));
    plot(xgrid, real(un), '-r','LineWidth',2);
    drawnow;
    %pause(0.01);
end




