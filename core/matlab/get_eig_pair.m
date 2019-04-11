clear all; curr_dir = cd; close all; clc;





gamma_vals = linspace(0,2/9,30); 
alpha_vals =  1:1:30;



p.gamma = gamma_vals(15);
p.alpha = alpha_vals(20);



file_name = ['gray_scott_evans_gamma_', ...
        num2str(round(p.gamma*1000000)), ...
        '_alpha_',num2str(round(p.alpha*1000000))];
ld = retrieve_it(curr_dir,'gray_scott',file_name,'data');
d = ld.var;
p = d.p;
lambda0 = d.root;

% controls
s.I = 30; % spatial infinity

% dependent controls
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

[s,e,m,c] = emcset(s,'front', LdimRdim(@A,s,p),'default');
% c.evans = @local_evans; % even Evans function
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


m.options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Refine',1,'Stats','off','MaxOrder',2);
m.ode_fun = @ode15s;

d = get_eig_fun(c,s,p,m,e,lambda0);

YRx = deval(d.Omega_R,0);
Omega = reshape(YRx(1:d.n*d.kr).',d.n,d.kr);
alpha_R = deval(d.alpha_R,0);
temp = Omega*alpha_R;
d.scl = temp(1);

x = linspace(s.L,-s.L,2001);
u = zeros(4,length(x));
for j = 1:length(x)
 u(:,j) = form_eig_func(x(j),d);
end

plot(x,u,'LineWidth',2);
h = xlabel('x');
set(h,'FontSize',22);
h = legend('u','w','v','z','Location','Best');
set(h,'FontSize',22);







