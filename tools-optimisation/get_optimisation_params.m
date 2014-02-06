function [const, deriv] = get_optimisation_params(f_trans, f, H)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 26/11/2013
% -------------------------------------------------------------------------
%
% GET_OPTIMISATION_PARAMS
% [const, deriv] = get_optimisation_params(f_trans, f, H)
%
% Given an energy function: E(f, H) = f_trans - H * f
% 
% get_optimisation_params outputs the coefficients of the energy function
% E(f, H) is linearised. We linearise by Taylor expansion:
%
%   E(f+df, H+dh) = E(f, H) + (dE/df')*df + (dE/dvecH')*dvecH
%
% This can be written in vector form as:
%
%   E(f+df, H+dh) = E(f, H) + [(dE/df') (dE/dvecH')] * [df; dvecH]
%
% Thus, the function get_optimisation_params returns the constant term
% (const = E(f, H)) and the derivative/Jacobian term (deriv = 
% [(dE/df') (dE/dvecH')]')
%
% Inputs:
%   - f_trans:  Position [x_trans; y_trans] of feature/point in transformed
%               co-ordinate frame B
%   - f:        Position [x; y] of feature/point in non-transformed
%               co-ordinate frame A
%   - H:        3x3 ransformation matrix from co-ordinate frame A to B, 
%               such that H * [f; 1] transforms feature point f to frame B.
%
% Outputs:
%   - const:    Constant term when E is linearised (2x1 vector)
%   - deriv:    Derivative term when E is linearised (10x2 matrix)

x = f(1);
y = f(2);

H_u = H(1:2,:);
H_l = H(3,:);

% Get constant term
const = f_trans - (H_u * [f; 1])/(H_l * [f; 1]);

% Get derivative/Jacobian term
deriv = [[ 1, 0, -x/(H(3,1)*x + H(3,2)*y + 1), 0, ...
    (x*(H(1,3) + H(1,1)*x + H(1,2)*y))/(H(3,1)*x + H(3,2)*y + 1)^2, ...
    -y/(H(3,1)*x + H(3,2)*y + 1), 0, ...
    (y*(H(1,3) + H(1,1)*x + H(1,2)*y))/(H(3,1)*x + H(3,2)*y + 1)^2, ...
    -1/(H(3,1)*x + H(3,2)*y + 1), 0]; ...
    ...
    [ 0, 1, 0, -x/(H(3,1)*x + H(3,2)*y + 1), ...
    (x*(H(2,3) + H(2,1)*x + H(2,2)*y))/(H(3,1)*x + H(3,2)*y + 1)^2, 0, ...
    -y/(H(3,1)*x + H(3,2)*y + 1), ...
    (y*(H(2,3) + H(2,1)*x + H(2,2)*y))/(H(3,1)*x + H(3,2)*y + 1)^2, 0, ...
    -1/(H(3,1)*x + H(3,2)*y + 1)]];

deriv = deriv';

end