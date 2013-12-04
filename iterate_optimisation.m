f_glob = [1; 1];
fs_loc{1} = [1; 0];
Hs{1} = eye(3);
fs_loc{2} = [1; -1];
Hs{2} = [cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1];
for i = 1:5
    [f_glob, Hs] = optimise_many_features(f_glob, fs_loc, Hs);
end

f_glob
Hs