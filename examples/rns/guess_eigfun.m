function out = guess_eigfun(x)
    out = [sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];
end

