function out = profile_eval(sol,x)

out = deval(sol,-abs(x));
out(2) = -sign(x)*out(2);
out(4) = -sign(x)*out(4);