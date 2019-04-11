function [eigfun, eigval] = eigenfunction()
    ld = load('eigfun.mat');
    eigfun = ld.p.eigfun;
    eigval = ld.p.eigfun.parameters;
end