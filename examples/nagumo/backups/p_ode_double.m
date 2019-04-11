function out = p_ode_double(x,y,P,n,lam)

temp = zeros(4,1);
for k = 0:n
    for l = 0:k
        if local_delta(n,k,l) == 1 
            temp1 = P{n-k+1}.fun(x);
            temp2 = P{k-l+1}.fun(x);
            temp3 = P{l+1}.fun(x);
            temp = temp + temp1.*temp2.*temp3;
        end
    end
end
qL = temp(1);
qR = temp(3);

BL = [0 1; 1+n*lam-6*sech(x)^2 0 ];
BR = [0 1; 1+n*lam-6*sech(-x)^2 0 ];

mu_L = sqrt(3*n+1);
mu_R = -sqrt(3*n+1);
Id = eye(2);

ind_L = 1:2;
ind_R = 3:4;

out = [(BL-mu_L*Id)*y(ind_L)+exp(-mu_L*x)*[0;-qL];
       (-BR+mu_R*Id)*y(ind_R)-exp(mu_R*x)*[0;-qR]];




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


