function [const, deriv] = get_optimisation_params_rep(f_trans, f, H)
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

xhat = f_trans(1);
yhat = f_trans(2);
H_tmp = num2cell(H(:));
[h11, h21, h31, h12, h22, h32, h13, h23, h33] = deal(H_tmp{:});
x = f(1);
y = f(2);

H_u = H(1:2,:);
H_l = H(3,:);

% Get constant term
const = f - (H_u * [f_trans; 1])/(H_l * [f_trans; 1]);

% Get derivative/Jacobian term
deriv = ...
    [[(h31*(h13 + h11*xhat + h12*yhat))/(h33 + h31*xhat + h32*yhat)^2 - h11/(h33 + h31*xhat + h32*yhat), ...
    (h32*(h13 + h11*xhat + h12*yhat))/(h33 + h31*xhat + h32*yhat)^2 - h12/(h33 + h31*xhat + h32*yhat), ...
    -xhat/(h33 + h31*xhat + h32*yhat), ...
    0, ...
    (xhat*(h13 + h11*xhat + h12*yhat))/(h33 + h31*xhat + h32*yhat)^2, ...
    -yhat/(h33 + h31*xhat + h32*yhat), ...
    0, ...
    (yhat*(h13 + h11*xhat + h12*yhat))/(h33 + h31*xhat + h32*yhat)^2, ...
    -1/(h33 + h31*xhat + h32*yhat), ...
    0]; ...
    ...
    [(h31*(h23 + h21*xhat + h22*yhat))/(h33 + h31*xhat + h32*yhat)^2 - h21/(h33 + h31*xhat + h32*yhat), ...
    (h32*(h23 + h21*xhat + h22*yhat))/(h33 + h31*xhat + h32*yhat)^2 - h22/(h33 + h31*xhat + h32*yhat), ...
    0, ...
    -xhat/(h33 + h31*xhat + h32*yhat), ...
    (xhat*(h23 + h21*xhat + h22*yhat))/(h33 + h31*xhat + h32*yhat)^2, ...
    0, ...
    -yhat/(h33 + h31*xhat + h32*yhat), ...
    (yhat*(h23 + h21*xhat + h22*yhat))/(h33 + h31*xhat + h32*yhat)^2, ...
    0, ...
    -1/(h33 + h31*xhat + h32*yhat)]];

deriv = deriv';

end