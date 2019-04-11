% Returns the right hand side of the equation y' = A*y + B
function out = ode_fun(x,y,P,n,lam,p)
    %0.075 each time. -> 0.017
    % Form A.
    tic;
    u0 = P{1}{1}(x);
    v0 = P{1}{3}(x);
    toc;
    tic;
    A = [0, 1, 0, 0; 
        1 - v0 ^ 2 - 3 * u0 ^ 2, 0, lam*n - 2 * u0 * v0 + 2*p.nu, 0;
        0, 0, 0, 1;
        lam*n - 2 * v0 * u0, 0, p.mu - 3 * v0^2 - u0^2, 0];
    
    toc;
    tic;
    % Form B, most of the slowdown is here.
    sm1 = 0;
    sm2 = 0;
    for j = 1:n
        Pj1 = P{j}{1}(x);
        Pj3 = P{j}{3}(x);
        
        for k = 1:n-(j-1)
            Pk1 = P{k}{1}(x);
            Pk3 = P{k}{3}(x);
            
            for l = 1:n-(j-1)-(k-1)
                if local_delta(n,j,k,l) == 1
                    Pl1 = P{l}{1}(x);
                    Pl3 = P{l}{3}(x);
              
                    sm1 = sm1 + Pj1*(Pk3*Pl3 + Pk1*Pl1);
                    %sm1 = sm1 + P{j}{1}(x) * P{k}{3}(x) * P{l}{3}(x);
                    %sm1 = sm1 + P{j}{1}(x) * P{k}{1}(x) * P{l}{1}(x);
                    
                    sm2 = sm2 + Pj3*(Pk3*Pl3 + Pk1*Pl1);
                    %sm2 = sm2 + P{j}{3}(x) * P{k}{3}(x) * P{l}{3}(x);
                    %sm2 = sm2 + P{j}{3}(x) * P{k}{1}(x) * P{l}{1}(x);
                end
            end
        end
    end   
    toc;
    tic;
    B = [0;sm1;0;sm2];
    
    % Form out.
    out = A*y + B;
    toc;
    fprintf("\n\n");

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


