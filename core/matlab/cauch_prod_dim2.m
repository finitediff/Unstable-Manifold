function out = cauch_prod_dim2(a,b)

% returns the coefficient of the Cauchy product of two vectors of size m by n

% out = sum(sum(a*flipud((flipud(b)).')));

out = sum(diag((fliplr(a)*fliplr(b.'))));



