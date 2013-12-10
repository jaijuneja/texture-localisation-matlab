function [f_glob_new, H_new] = optimise_single_feature(f_glob, f_loc, H)

[a, b] = get_optimisation_params(f_glob, f_loc, H);

% Calculate parameters for quadprog
% Recall the formula:
% Energy = a'*a + 2*a'*b'*delta + delta'*b*b'*delta
% Quadprog takes form quadprog(G, h) for 1/2*x'*G*x + h'*x

E_before = a' * a

G = 2 * (b * b');
h = 2 * a' * b';
h = h';

delta = quadprog(G, h);

f_glob_new = f_glob + delta(1:2);

H_new = H;
for i = 1:8
    H_new(i) = H(i) + delta(2+i);
end

E_after = a'*a + 2*a'*b'*delta + delta'*(b*b')*delta

end