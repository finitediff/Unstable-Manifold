function out = local_evans(yl,yr,lambda,s,p,m,e)

fun=str2func(e.evans);
out=fun(yl,yr,lambda,s,p,m,e);

%--------------------------------------------------------------------------
% reg_reg_polar
%--------------------------------------------------------------------------

function out = reg_reg_polar(WL,WR,lambda,s,p,m,e)

OmegaL0 = orth(WL);
alphaL = OmegaL0'*WL;
muL = trace(OmegaL0'*e.LA(e.Li(1),lambda,s,p)*OmegaL0);
[omegal,gammal] = manifold_polar(e.Li,OmegaL0,lambda,e.LA,s,p,m,e.kl,muL);

E = [1 0; 
     0 0;
     0 1;
     0 0];

out = (det(alphaL)*gammal)*det([omegal,E]);



