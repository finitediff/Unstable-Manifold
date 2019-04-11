
function [out,s,c,m,e] = evan_root_solver2(p,s,tol)


d.p = p;

% dependent controls
s.L = -s.I;
s.R = s.I;

% solve the profile
ode_fun = @(x,y,params)(profile_ode(x,y,p));
bc_fun = @(ya,yb,params)(profile_bc(ya,yb,p));
x = linspace(s.L,0,30);
guess = @(x)(profile_guess(x,p)); % NOT USING CONTINUATION
solinit = bvpinit(x,guess);
options = bvpset('RelTol',1e-12,'AbsTol',1e-12);
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

% plot the profile
% plot_profile(p,s);

R = high_frequency_bound(p,s);
% fprintf('R: %4.4g\n',R);
d.R = R;

d.s = s;


% 
% Evans function
% 


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


m.options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Refine',1,'Stats','off','MaxOrder',2);
m.ode_fun = @ode15s;

% [w,dom] = contour(c,s,p,m,e,linspace(1e-2,R,30))
% plot(dom,w,'.-k');


fun = @(a,b)(contour(c,s,p,m,e,linspace(a,b,3)));

out = evan_root_bisector(fun,1e-2,R,tol);





















