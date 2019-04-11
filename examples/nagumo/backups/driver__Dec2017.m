clc; clear all; close all; beep off;

% Nagumo - reaction diffusion equation
% 
% u_t = u_{xx} - u + u^3

X = 30;
n_vals = 2:1:3;


XL = -X;
XR = X;

% parametrization
p.c = 1;

P{1} = @(x)(sqrt(2)*sech(x));
P{2} = @(x)(eig_fun(x));

lam = 3; % eigenvalue of the linearized PDE


cnt = 0;

for n = n_vals
    
    n
    

    ode_handle = @(x,y)(p_ode(x,y,P,n,lam));

    bc_handle = @(ya,yb)(bc(ya,yb,n));

    options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000);

    pre_guess = @(x)([sech(x);-sech(x)*tanh(x)]);

    solinit = bvpinit([XL,XR],pre_guess);

    sol = bvp5c(ode_handle,bc_handle,solinit);
    
    cnt = cnt + 1;
    d{cnt} = sol;
    
    P{n+1} = @(x)(interp(x,sol));

end


figure;
hold on;
x = linspace(XL,XR,101);
for j = 3:n_vals(end)
    y = zeros(size(x));
    for k = 1:length(x)
        y(k) = P{j}(x(k));
    end
    plot(x,y,'-k');
    pause(1)
end

