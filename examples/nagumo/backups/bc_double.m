function out = bc_double(ya,yb,n)


PL = [-sqrt(1+3*n),1,0,0];
PR = [0,0,sqrt(1+3*n),1];

out = [PL*ya;PR*ya;yb(1:2)-yb(3:4)];











