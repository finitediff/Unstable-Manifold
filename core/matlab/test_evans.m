
s.I = 15;
p.alpha = 6;
p.gamma = 1/9;

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
c.evans = @local_evans; % even Evans function
c.basisL = @analytic_basis_local;
c.basisR = @analytic_basis_local;
% c.single_process = 1;

% 
% process and display data
%

% display a waitbar
c.stats = 'off';
c.refine = 'on';
c.ksteps = 2^8;
c.lambda_steps = 0;
c.tol = 0.2;


m.options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Refine',1,'Stats','off','MaxOrder',2);
m.ode_fun = @ode15s;


%
% preimage contour
%

% circpnts=30; imagpnts=30; spread=2; inner_pnts = 10; inner_radius = 10^(-2);
% preimage = semicirc_left_inclusion(circpnts,imagpnts,inner_pnts,c.ksteps,R, ...
%     spread,inner_radius,c.lambda_steps);

circpnts=30; imagpnts=30; spread=2; inner_pnts = 10; inner_radius = 0.01;
preimage=semicirc2(circpnts,imagpnts,inner_pnts,c.ksteps,R, ...
    spread,inner_radius,c.lambda_steps);
d.inner_radius = inner_radius;

%
% compute Evans function
%

tic
halfw = contour(c,s,p,m,e,preimage);
d.time1 = toc;
halfw = halfw/halfw(1);
w = [halfw fliplr(conj(halfw))];
d.wnd1 = winding_number(w);
d.halfw1 = halfw;


preimage = half_circle(circpnts,c.ksteps,inner_radius,c.lambda_steps);
tic
halfw = contour(c,s,p,m,e,preimage);
d.time2 = toc;
halfw = halfw/halfw(1);
w = [halfw fliplr(conj(halfw))];
d.wnd2 = winding_number(w);
d.halfw2 = halfw;






















