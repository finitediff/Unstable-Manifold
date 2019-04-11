function [out, projects] = analytic_basis_local(projection,x,preimage,s,p,A,posneg,eps,Q,p_old_in)
% [out, projects] = analytic_basis(projection,x,preimage,s,p,A,posneg,eps,Q,p_old_in)
%
% Returns an analytic basis at specified infinity using the method of Kato.
%
% Input "projection" is a function handle to the projection function to be
% used, "x" is the numerical value of infinity, "preimage" is the contour
% on which the Evans function is computed, "s" and "p" are structures
% explained in the STABLAB documentation, "A" is a function handle to the
% Evans matrix, "posneg" is 1 or -1 determining which space the
% projection function should return, and "eps" is the tolerance in the
% projection function. If input Q and p_old_in are specified, then the
% first value of the analtyic basis is set to Q with projection p_old_in.
% This allows for continuing or filling in a previous analtyic basis
% computation. When these two inputs are not specified, they are chosen
% automtically.

iterations = size(preimage,2);
[p_old, Q1] = projection(A(x,preimage(1),s,p),posneg,eps);
[n,k]=size(Q1);
out = zeros(n,k,iterations);
projects = zeros(size(p_old,1),size(p_old,2),iterations);

for j=1:iterations
    
    mu1 = -sign(x)*sqrt(preimage(j)+p.alpha);
    mu2 = -sign(x)*sqrt(preimage(j)+1/p.gamma);
    
%     out(:,1,j) = [1;mu1;0;0]/sqrt(1+mu1^2);
%     out(:,2,j) = [0;0;1;mu2]/sqrt(1+mu2^2);

    out(:,1,j) = [1;mu1;0;0];
    out(:,2,j) = [0;0;1;mu2];

end
 



