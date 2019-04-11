function out = cauchy_prod3_smaller_terms(a,b,c)

a = squeeze(a);
b = squeeze(b);
c = squeeze(c);

[m1,n1] = size(a);
m = m1-1;
n = n1-1;


out = 0;


for ind = 0:m
    for j = 0:ind
        for k = 0:n
            for l = 0:k
                
                if (ind==0)&&(k==0)
                   continue 
                end
                if (ind==m)&&(k==n)&&(j==0)&&(l==0)
                   continue 
                end
                if (j==m)&&(l==n)
                    continue
                end
                
                out = out + a(m-ind+1,n-k+1)*b(ind-j+1,k-l+1)*c(j+1,l+1);
                
            end
        end
    end
end




