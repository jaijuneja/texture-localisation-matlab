function [const, deriv] = constrain_scale(H)

A = H(1:2,1:2);

const = det(A) - 1;

ddetA_dA = det(A) * inv(A)';

deriv = 2 * const * ddetA_dA(:);

const = const^2;

end