function out = profile_bc(ya,yb,p)

mu1 = -sqrt(p.alpha);
mu2 = -sqrt(p.gamma);

ya0 = [1;0;0;0];

vec1 = [mu1,1,0,0];
vec2 = [0,0,mu2,1];

out = [vec1*(ya(1:4)-ya0);
    vec2*(ya(1:4)-ya0);
    yb(2);
    yb(4)];












