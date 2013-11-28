function [world, cor] = optimise_features(world, cor)
% Stuff goes here
end
%%
syms dh11 dh12 dh13 dh21 dh22 dh23 dh31 dh32 dh33 x y;
dH = [
     dh11 dh12 dh13
     dh21 dh22 dh23
     dh31 dh32 dh33
     ];
dH_u = dH(1:2,:);
dH_l = dH(3,:);
f_loc = [x; y; 1];

f_glob = @(dH_u, dH_l, f_loc)(dH_u * f_loc)/(dH_l * f_loc);

f = f_glob(dH_u, dH_l, f_loc);

dvecH = sym('dvecH', [9 1]);
dx_dvecH = sym('dx_dvecH', [1 9]);
dy_dvecH = sym('dy_dvecH', [1 9]);

for i = 1:9
    dvecH(i) = dH(i);
    dx_dvecH(i) = diff(f(1), dH(i));
    dy_dvecH(i) = diff(f(2), dH(i));
end

df_dH = [dx_dvecH; dy_dvecH];

% After evaluating df_dH, f_star and f, do this:
% A = 2 * (df_dH' * df_dH);
% syms x_star y_star;
% f_star = [x_star; y_star];
% b = 2 * (df_dH' * (f_star - f)');
% quadprog(A, b);