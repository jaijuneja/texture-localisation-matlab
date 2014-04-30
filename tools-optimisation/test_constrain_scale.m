A = rand(2,2);
startval = -1;
stepsize = 0.1;
endval = 1;

iters = (endval - startval)/stepsize + 1;
truevals = zeros(1, iters);
estvals = zeros(1, iters);
iter = 0;
vals = startval : stepsize : endval;
for step = vals
    iter = iter + 1;
    dA = repmat(step, 2, 2);
    truevals(iter) = (det(A+dA)-1)^2;
    
    [const, deriv] = constrain_scale(A);
    estvals(iter) = const + deriv * dA(:);
end

plot(vals, truevals, 'r', vals, estvals, 'b');