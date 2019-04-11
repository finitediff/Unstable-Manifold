function out = p_ode(x,y,P,n,lam)

sm = 0;
for k = 0:n
    for l = 0:k
        if local_delta(n,k,l) == 1            
            sm = sm + P{n-k+1}(x)*P{k-l+1}(x)*P{l+1}(x);
        end
    end
end

B = [0 1; 1+n*lam-6*sech(x)^2 0 ];

if x <= 0
   mu = sqrt(3*n+1)*eye(2); 
elseif x > 0
   mu = sqrt(3*n+1)*eye(2);
end

out = (B-mu*eye(2))*y+exp(-mu*x)*[0; -sm];



% delta function
function out = local_delta(n,k,l)


if (k == 0) && (l == 0)
    out = 0;
elseif (k == n) && (l == n)
    out = 0;
elseif (k == n) && (l == 0)
    out = 0;
else
   out = 1;
end


