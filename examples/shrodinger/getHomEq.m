function P = getHomEq(x, P, n, L, ode_fun, bc, lam, p, bvpPoints)
    %% Find the homological equations
    options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);

    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam,p));
    bc_handle = @(ya,yb)(bc(ya,yb,n,lam,p)); 
    total_fun = @(x) P{n}(x);
    
    solinit = bvpinit(linspace(-L,L,bvpPoints),total_fun);
    sol = bvp5c(ode_handle,bc_handle,solinit,options);  
    temp = deval(sol,x);

    % Convert into chebychev polynomial.
    for j = 1:4 
        [cf, fun{j}] = get_chebyshev_coefficients(temp(j,:),-L,L,'first kind');
    end
    P{n+1} = @(x) [fun{1}(x).'; fun{2}(x).'; fun{3}(x).'; fun{4}(x).'];
    
end

