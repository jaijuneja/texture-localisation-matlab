function [f_glob_new, H_new] = optimise_many_features(f_glob, fs_loc, Hs)
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
    [a, b] = get_optimisation_params(f_glob, fs_loc{i}, Hs{i});
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
H_new = Hs;
for i = 1:numImagesInMap
    imageNdx = numGlobalFeats*2+i*8-7;
    for j = 1:8
        H_new{i}(j) = H_new{i}(j) + delta(imageNdx+j-1);
    end
end

E_after = A'*A + h'*delta + delta'*G*delta

% f_glob = [1; 1];
% fs_loc{1} = [1; 0];
% Hs{1} = eye(3);
% fs_loc{2} = [1; -1];
% Hs{2} = [cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1];

figure
p1 = subplot(1,2,1);
plot(p1, f_glob(1), f_glob(2), 'rx'); hold on

for j = 1:numLocalFeats
    locfeat = Hs{j} * [fs_loc{j}; 1];
    locfeat = locfeat(1:2)/locfeat(3);
    plot(p1, locfeat(1), locfeat(2), 'bx');
end
hold off

p2 = subplot(1,2,2);
plot(p2, f_glob_new(1), f_glob_new(2), 'rx'); hold on
for j = 1:numLocalFeats
    locfeat_new = H_new{j} * [fs_loc{j}; 1];
    locfeat_new = locfeat_new(1:2)/locfeat_new(3);
    plot(p2, locfeat_new(1), locfeat_new(2), 'bx');
end
hold off

end