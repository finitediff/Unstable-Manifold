function y0 = eval_P(P,M,N,sig1,sig2)

y0 = zeros(4,1);
for m = 0:M
    for n = 0:N
        y0 = y0+P(:,m+1,n+1)*sig1^m*sig2^n;
    end
end