function y = form_eig_func(x,d)
% system for non-conservative form

% evaluate solution to problem
if x >=0
    YRx = deval(d.Omega_R,x);
    Omega = reshape(YRx(1:d.n*d.kr).',d.n,d.kr);
    alpha_R = deval(d.alpha_R,x);
    y = Omega*alpha_R;
else
    YLx = deval(d.Omega_L,x);
    Omega = reshape(YLx(1:d.n*d.kl).',d.n,d.kl);
    alpha_L = deval(d.alpha_L,x);
    y = Omega*alpha_L;
end

y = y/d.scl;
