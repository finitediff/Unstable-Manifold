function out = profile_ode(x,y,p)

% RHS of profile ODE equations

u = y(1);
w = y(2);
v = y(3);
z = y(4);

out = [ w; (u*v^2-p.alpha*(1-u)); z; (v-u*v^2)/p.gamma];











