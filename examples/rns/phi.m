function [exp_A, exp_A_prime] = phi(p,e_n)
    %Load parameters.
    cnu = p.cnu;
    Ti = p.Ti_weight;
    
    exp_A = zeros(size(e_n));
    for j = 1:length(e_n)
       if cnu*e_n(j)-Ti > 0
           exp_A(j) = exp(-A./(-Ti+e_n(j)/cnu));
           exp_A_prime(j) = 1/cnu*exp(A./(-Ti+e_n(j)/cnu)^2);
       else
           exp_A(j) = 0;
           exp_A_prime(j) = 0;
       end
    end
end

