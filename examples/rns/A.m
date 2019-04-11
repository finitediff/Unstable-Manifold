function G = A(x,lambda,s,p)
% G = A(x,lambda,s,p) 
%
% Regular Evans matrix for rNS

% profile
%temp = soln(x,s);
x
temp = s.profile_fun(x);
u = temp(1);
e = temp(2);
z = temp(3);
tau = p.tau_minus - (1/p.s)*(u - p.u_minus);

% profile derivatives
%der = profile_ode(0,temp,s,p);
der = s.profile_der(x, 10e-10);
u_x = der(1);
e_x = der(2);
z_x = der(3);

% phi and phi'
[phi_e, Dphi_e] = phi(p,e);

% partials of pressure function
p_tau = -p.g *e/tau^2;
p_e = p.g/tau;
p_z = 0; 
p_p = p.g*e/tau;

% entries of matrix G
g11 = lambda/p.s;
g12 = [0 0 0];
g13 = [lambda/p.s 0 0];

g21 = zeros(3,1);
g22 = zeros(3,3);
g23 = [ 
                -lambda           0                                               0;
                -lambda*u       -lambda+p.q*Dphi_e*z/p.c       p.q*phi_e;
                0                     -Dphi_e*z/p.c                           -lambda-phi_e
           ];
       
  g31 = (-1/p.s)*[ 
                            tau*p_tau/p.nu+u_x/tau;
                            e_x/tau; 
                            2*z_x/tau
                         ];
  g32 = - [ 
                    tau/p.nu                    0                               0
                    -tau*u*p.c/p.kappa    tau*p.c/p.kappa       0
                    0                               0                               tau^2/p.d
                ];
            
C1 = (-1/p.s)*(p_tau+p.nu*u_x/tau^2)-p.s;
C2 = (-1/p.s)*(p_tau*u+p.nu*u*u_x/tau^2 +p.kappa*e_x/(tau^2*p.c))-p.s*u+p_p-p.nu*u_x/tau;
C3 = (-1/p.s)*(2*p.d*z_x/tau^3);
            
g33 = [  
                tau*C1/p.nu                                 tau*p_e/p.nu                                                        tau*p_z/p.nu;
                -tau*u*C1*p.c/p.kappa+tau*C2*p.c/p.kappa    -p_e*tau*u*p.c/p.kappa-p.s*tau*p.c/p.kappa+tau*u*p_e*p.c/p.kappa    0;
                tau^2*C3/p.d                                0                                                                   -p.s*tau^2/p.d
            ];
        
% Evans matrix G
G = [ 
            g11 g12 g13;
            g21 g22 g23;
            g31 g32 g33
        ]; % 7 x 7?
    size(G)
% G = [ 
%             0   1 0   0 0   0;
%             g11 0 g12 0 g13 0;
%             0   0 0   1 0   0;
%             g21 0 g22 0 g23 0;
%             0   0 0   0 0   1;
%             g31 0 g32 0 g33 0
%        ];

