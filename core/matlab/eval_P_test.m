function y0 = eval_P_test(P,M,N,sig1,sig2,mu1_u,mu2_u)

y0 = zeros(4,1);
for m = 0:M
    for n = 0:N
        y0 = y0+(m*mu1_u+n*mu2_u)*P(:,m+1,n+1)*sig1^m*sig2^n;
    end
end