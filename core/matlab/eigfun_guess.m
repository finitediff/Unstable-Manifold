function u = eigfun_guess(x,d)



u = zeros(4,length(x));
for j = 1:length(x)
 u(:,j) = form_eig_func(x(j),d)/d.u_scl;
end







