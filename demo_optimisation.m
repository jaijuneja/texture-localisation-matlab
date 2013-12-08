%% Simple single feature, multiple view scenario

f_glob = [1; 1];
f_loc{1} = [1; 0];
H{1} = eye(3);
f_loc{2} = [1; -1];
H{2} = [cos(pi/3) -sin(pi/3) 0; sin(pi/3) cos(pi/3) 0; 0 0 1];

numLocalFeats = length(f_loc);
numIterations = 5;
locfeats = cell(1, numLocalFeats);

for i = 1:numIterations
    subplot(1,numIterations,i)
    plot(f_glob(1), f_glob(2), 'rx'); hold on
    [f_glob, H] = optimise_many_features(f_glob, f_loc, H);
    
    for j = 1:numLocalFeats
    locfeat = H{j} * [f_loc{j}; 1];
    locfeat = locfeat(1:2)/locfeat(3);
    plot(locfeat(1), locfeat(2), 'bx');
        if isequal(i, numIterations)
            locfeats{j} = locfeat;
        end
    end
end

%% Real application
clear
load('/Users/jai/Documents/MATLAB/4yp/test_images/wall/etc/index.mat')
load('/Users/jai/Documents/MATLAB/4yp/test_images/wall/etc/cor.mat')
load('/Users/jai/Documents/MATLAB/4yp/test_images/wall/etc/world.mat')

%% Make a sample world and cor struct
% Some basic cases first:
% e.g. only move global features, only adjust H etc.
%% Global feature bundle adjustment
numIters = 5;
for i = 1:numIters
    [world, cor] = bundle_adjustment(world, cor, 'perspDistPenalty', 1e5, ...
        'onlyOptimiseH', true);
end
%% Local feature bundle adjustment
numIters = 5;
for i = 1:numIters
    [world, cor] = bundle_adjustment_local(world, cor, 'perspDistPenalty', 1e5, ...
        'onlyOptimiseH', false);
end
%%
tic
world = build_world(index, cor);
toc
%%
plot_features(world)

%%
% figure;
hold on
% plot_feature_matches(world, cor);
plot_features(world);
hold on
plot_transformations(cor, index);