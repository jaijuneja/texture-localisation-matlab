function [a, b] = get_optimisation_params(f_glob, f_loc, H)
x = f_loc(1);
y = f_loc(2);

H_u = H(1:2,:);
H_l = H(3,:);

a = f_glob - (H_u * [f_loc; 1])/(H_l * [f_loc; 1]);

b = [[ 1, 0, -x/(H(3,1)*x + H(3,2)*y + 1), 0, ...
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

b = b';

end

%% Calculation of derivative terms

% syms h11 h12 h13 h21 h22 h23 h31 h32 h33 x y x_star y_star;
% H = [
%      h11 h12 h13
%      h21 h22 h23
%      h31 h32 1
%      ];
% H_u = H(1:2,:);
% H_l = H(3,:);
% f_loc = [x; y; 1];
% f_star = [x_star; y_star];
% 
% getOmega = @(H_u, H_l, f_loc, f_star)(f_star - (H_u * f_loc)/(H_l * f_loc));
% 
% omega = getOmega(H_u, H_l, f_loc, f_star);
% 
% dvecH = sym('dvecH', [8 1]);
% dx_dvecH = sym('dx_dvecH', [1 8]);
% dy_dvecH = sym('dy_dvecH', [1 8]);
% syms dx_dxstar dx_dystar dy_dxstar dy_dystar;
% 
% for i = 1:8
%     dvecH(i) = H(i);
%     dx_dvecH(i) = diff(omega(1), H(i));
%     dy_dvecH(i) = diff(omega(2), H(i));
% end
% 
% dx_dxstar = diff(omega(1), x_star);
% dx_dystar = diff(omega(1), y_star);
% dy_dxstar = diff(omega(2), x_star);
% dy_dystar = diff(omega(2), y_star);
% 
% b = [dx_dxstar dx_dystar dx_dvecH; dy_dxstar dy_dystar dy_dvecH];