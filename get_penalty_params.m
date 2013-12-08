function [a31, b31, a32, b32] = get_penalty_params(H, penalty)

% penalty*(h31 + dh31)^2 = penalty * (h31^2 + 2h31*dh31 + dh31^2)

a31 = H(3,1) * sqrt(penalty);
a32 = H(3,2) * sqrt(penalty);
b31 = 2 * H(3,1) *  penalty;
b32 = 2 * H(3,2) *  penalty;

end