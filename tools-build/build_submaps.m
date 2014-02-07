function submaps = build_submaps(world, varargin)

opts.numSubmaps = 10;
opts.numRepetitions = 5;
opts.featsPerSubmap = false;
opts = vl_argparse(opts, varargin);

feats_glob = world.features_global(3:4, :);
num_feats = length(feats_glob);

if opts.featsPerSubmap
    opts.numSubmaps = ceil(num_feats / opts.featsPerSubmap);
end

% Get submap centres and indexes of features belonging to each submap using
% k-means clustering
[sm_ctrs, sm_idxs] = vl_kmeans(feats_glob, opts.numSubmaps, ...
    'NumRepetitions', opts.numRepetitions);

% Create matrix of colours
ColOrd = get(gca,'ColorOrder');
[m,~] = size(ColOrd);

initvals = cell(1, opts.numSubmaps);
submaps = struct(...
    'center', initvals, ...
    'H_to_world', initvals, ...
    'feats', initvals, ...
    'words', initvals);

for i = 1:opts.numSubmaps
    ColRow = rem(i,m);
    if ColRow == 0
      ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);
    
    submaps(i).center = sm_ctrs(:, i);
    submaps(i).H_to_world = [1 0 sm_ctrs(1, i); 0 1 sm_ctrs(2, i); 0 0 1];
    feats = world.features_global(:, sm_idxs==i);
    feats(3, :) = feats(3, :) - sm_ctrs(1, i);
    feats(4, :) = feats(4, :) - sm_ctrs(2, i);
    submaps(i).feats = feats;
    submaps(i).words = world.words_global(:, sm_idxs==i);
    
    plot(feats_glob(1, sm_idxs==i), feats_glob(2, sm_idxs==i), '.', ...
        'Color', Col, 'MarkerSize',12)
    hold on
end
plot(sm_ctrs(1,:),sm_ctrs(2,:),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(sm_ctrs(1,:),sm_ctrs(2,:),'ko',...
     'MarkerSize',12,'LineWidth',2)
set(gca, 'YDir', 'reverse')
axis equal, axis tight
end