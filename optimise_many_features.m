function [f_glob_new, H_new] = optimise_many_features(f_glob, f_loc, H)
% Test:
% f_glob = [1; 1];
% fs_loc{1} = [1; 0];
% Hs{1} = eye(3);
% fs_loc{2} = [1; -1];
% Hs{2} = [cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1];

% Pre-populate global parameters
% Need to use sparse matrices for the real deal!
numGlobalFeats = 1;
numLocalFeats = 2;
numImagesInMap = 2;
A = zeros(numLocalFeats * 2, 1);
h = zeros(numGlobalFeats * 2 + numImagesInMap * 8, 1);
G = zeros(numGlobalFeats * 2 + numImagesInMap * 8);
imageID = 0;
for i = 1:numLocalFeats
    globalFeatID = 1;
    featNdx = 2*globalFeatID-1;
    imageID = imageID + 1; % Need an if statement to check image ID
    imageNdx = numGlobalFeats*2+imageID*8-7;
    [a, b] = get_optimisation_params(f_glob, f_loc{i}, H{i});
    A(2*i-1:2*i) = a;
    h_Loc = 2 * a' * b';
    h(featNdx:featNdx+1) = h(featNdx:featNdx+1) + h_Loc(1:2)';
    h(imageNdx:imageNdx+7) = h(imageNdx:imageNdx+7) + h_Loc(3:10)';
    G_Loc = (b * b');
    G(featNdx:featNdx+1,featNdx:featNdx+1) = ...
        G(featNdx:featNdx+1,featNdx:featNdx+1) + G_Loc(1:2,1:2);
    G(featNdx:featNdx+1,imageNdx:imageNdx+7) = ...
        G(featNdx:featNdx+1,imageNdx:imageNdx+7) + G_Loc(1:2,3:10);
    G(imageNdx:imageNdx+7,featNdx:featNdx+1) = ...
        G(imageNdx:imageNdx+7,featNdx:featNdx+1) + G_Loc(3:10,1:2);
    G(imageNdx:imageNdx+7,imageNdx:imageNdx+7) = ...
        G(imageNdx:imageNdx+7,imageNdx:imageNdx+7) + G_Loc(3:10,3:10);
end

% Calculate parameters for quadprog
% Recall the formula:
% Energy = a'*a + 2*a'*b'*delta + delta'*b*b'*delta
% Quadprog takes form quadprog(G, h) for 1/2*x'*G*x + h'*x

E_before = A' * A

G = 2 * G;

delta = quadprog(G, h);

ndx = 1;
f_glob_new = f_glob + delta(ndx:ndx+1);
H_new = H;
for i = 1:numImagesInMap
    imageNdx = numGlobalFeats*2+i*8-7;
    for j = 1:8
        H_new{i}(j) = H_new{i}(j) + delta(imageNdx+j-1);
    end
end

%delta(3:end) = 0 ;

E_after = A'*A + h'*delta + 0.5*delta'*G*delta;

if 1
current = 0;
for i = 1:numLocalFeats
    img_ndx = 2 + 8 * i - 7;
    [a, b] = get_optimisation_params(f_glob, f_loc{i}, H{i});
    delta_new = [delta(1:2); delta(img_ndx:img_ndx+7)];
    current = current + norm(a + b' * delta_new)^2;
end
assert(abs(current - E_after) < 1e-5) ;
end

end