clear all; close all; clc; beep off; curr_dir = cd;

%
% parameters
%

p.gamma = 1/9;
p.alpha = 1/p.gamma;

% controls
s.I = 30; % spatial infinity


tol = 1e-3;
[lambda0,s,c,m,e] = evan_root_solver(p,s,tol);

d = get_eig_fun(c,s,p,m,e,lambda0);

YRx = deval(d.Omega_R,0);
Omega = reshape(YRx(1:d.n*d.kr).',d.n,d.kr);
alpha_R = deval(d.alpha_R,0);
temp = Omega*alpha_R;
d.scl = temp(1);

x = linspace(s.L,-s.L,2001);
u = zeros(4,length(x));
u_scl = norm(form_eig_func(0,d));
d.u_scl = u_scl;
for j = 1:length(x)
 u(:,j) = form_eig_func(x(j),d)/u_scl;
end

hold on;
plot(x,u,'LineWidth',2);
h = xlabel('x');
set(h,'FontSize',22);
h = legend('u','w','v','z','Location','Best');
set(h,'FontSize',22);


guess = @(x)(eigfun_guess(x,d));
ode_fun = @(x,y,lambda)(A(x,lambda,s,p)*y);
fun_bc = @(ya,yb,lambda)(bc_fun(ya,yb,lambda,p));

solinit = bvpinit([s.L,0],guess,lambda0);
options = bvpset('RelTol',1e-10,'AbsTol',1e-10);
eigfun_sol = bvp5c(ode_fun,fun_bc,solinit,options);

plot(eigfun_sol.x,eigfun_sol.y,'--g')


data.p = p;
data.profile_sol = s.sol;
data.eigfun_sol = eigfun_sol;
data.lambda0 = eigfun_sol.parameters;

file_name = ['gray_scott_eigfun_gamma_', ...
             num2str(round(p.gamma*1000000)), ...
             '_alpha_',num2str(round(p.alpha*1000000))];
saveit(curr_dir,'gray_scott',data,file_name,'eigfun');

return




 
% 
% process and display data
%

d = evans_batch(p,s);

w = [d.halfw fliplr(conj(d.halfw))];
fprintf('Winding Number: %1.1d\n',d.wnd);
plot_evans(w);
