function plot_profile(p,s)

% plot the profile

x = linspace(s.L,0,200);
u = zeros(size(x));
w = u;
v = u;
z = u;

for j = 1:length(x)
    temp = deval(s.sol,x(j));
    u(j) = temp(1);
    w(j) = temp(2);
    v(j) = temp(3);
    z(j) = temp(4);
end

figure;
hold on;
plot(x,u,'-k','LineWidth',2);
plot(x,w,'-.k','LineWidth',2);
plot(x,v,'--k','LineWidth',2);
plot(x,z,':k','LineWidth',2);

h = legend('u','w','v','z','Location','Best');
set(h,'FontSize',22);
h = xlabel('x');
set(h,'FontSize',22);
h = ylabel('profile');
set(h,'FontSize',22);
h = gca;
set(h,'FontSize',22);
