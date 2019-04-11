% Returns the right hand side of the equation y' = A*y + B
function out = ode_fun(x,y,P,n,lam,p)

    %0.075 each time. -> 0.017 -> 0.002
    
    % Deval all the terms beforehand to drastically speed it up.
    P1 = zeros(1,n);
    P3 = zeros(1,n);
    for j = 1:n
        temp = P{j}(x);
        P1(j) = temp(1);
        P3(j) = temp(3);
    end
    
    % Form A.
    u0 = P1(1);
    v0 = P3(1);
    A = [0, 1, 0, 0; 
        1 - v0 ^ 2 - 3 * u0 ^ 2, 0, lam*n - 2 * u0 * v0 + 2*p.nu, 0;
        0, 0, 0, 1;
        -lam*n - 2 * v0 * u0, 0, p.mu - 3 * v0^2 - u0^2, 0];
    
    % Form B, most of the slowdown is here.
    sm1 = 0;
    sm2 = 0;
    
    % j-1+k-1+l-1 = n
    for j = 0:n       
        for k = 0:n-j
            l = n-j-k;
            if local_delta(n,j,k,l) == 1
                secondTerm = P3(k+1)*P3(l+1) + P1(k+1)*P1(l+1);
                sm1 = sm1 + P1(j+1)*secondTerm;
                sm2 = sm2 + P3(j+1)*secondTerm;
            end
        end
    end   
    B = [0;sm1;0;sm2];
    
    % Form out.
    out = A*y - B;
    %fprintf("\n\n");

% delta function, returns 1 if j,k,l all not equal to n.
function out = local_delta(n,j,k,l)

    if (j == n)
        out = 0;
    elseif (k == n)
        out = 0;
    elseif (l == n)
        out = 0;
    else
       out = 1;
    end


