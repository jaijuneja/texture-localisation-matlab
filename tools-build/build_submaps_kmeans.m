function submaps_new = build_submaps_kmeans(model, world, varargin)

opts.numCenters = 15;
opts.numRepetitions = 5;
opts.featsPerSubmap = [];
opts.rerankDepth = 5;
opts.showPlot = true;
opts = vl_argparse(opts, varargin);

feats_glob = world.features_global(3:end, :);
num_feats = length(feats_glob);

if ~isempty(opts.featsPerSubmap)
    opts.numCenters = ceil(3 * num_feats / opts.featsPerSubmap);
end

% Get submap centres and indexes of features belonging to each submap using
% k-means clustering
[sm_ctrs, sm_idxs] = vl_kmeans(feats_glob(1:2,:), opts.numCenters, ...
    'NumRepetitions', opts.numRepetitions, 'Initialization', 'PLUSPLUS');

% Create matrix of colours
ColOrd = get(gca,'ColorOrder');
[m,~] = size(ColOrd);

initvals = cell(1, opts.numCenters);

submaps.center = initvals;
submaps.H_to_world = initvals;
submaps.feats = initvals;
submaps.words = initvals;
submaps.ids = 1:opts.numCenters;

for i = 1:opts.numCenters
    ColRow = rem(i,m);
    if ColRow == 0
      ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);
    
    submaps.center{i} = sm_ctrs(:, i);
    submaps.H_to_world{i} = [1 0 sm_ctrs(1, i); 0 1 sm_ctrs(2, i); 0 0 1];
    feats = world.features_global(3:end, sm_idxs==i);
    feats(1, :) = feats(1, :) - sm_ctrs(1, i);
    feats(2, :) = feats(2, :) - sm_ctrs(2, i);
    submaps.feats{i} = feats;
    submaps.words{i} = world.words_global(:, sm_idxs==i);
    
    if opts.showPlot
        plot(feats_glob(1, sm_idxs==i), feats_glob(2, sm_idxs==i), '.', ...
            'Color', Col, 'MarkerSize',12)
        hold on
    end
end

% Compute Delaunay triangulation
ctrs = cell2mat(submaps.center)';
DT = delaunayTriangulation(ctrs);

% Compute map bounding region
xmax = max(world.features_global(3,:));
xmin = min(world.features_global(3,:));
ymax = max(world.features_global(4,:));
ymin = min(world.features_global(4,:));
region = [xmin xmax ymin ymax];

% Re-calculate Delaunay triangulation given the clipped region
[DT, numNewCtrs] = clip_voronoi(DT, region);
DT = struct('Points', DT.Points, 'ConnectivityList', DT.ConnectivityList, ...
    'Constraints', DT.Constraints);

% Remove the new centers that were added for clipping purposes
DT.Points(1:numNewCtrs, :) = [];
toRemove = ismember(DT.ConnectivityList, 1:numNewCtrs);
[toRemove, ~] = find(toRemove);
toRemove = unique(toRemove);
DT.ConnectivityList(toRemove, :) = [];
DT.ConnectivityList = DT.ConnectivityList - numNewCtrs;

% Expand submaps to include Delaunay neighbours
numSubmaps = size(DT.ConnectivityList, 1);
numWords = size(model.index.histograms, 1);
histograms = zeros(numWords, numSubmaps);

initvals = cell(1, numSubmaps);

submaps_new.center = initvals;
submaps_new.H_to_world = initvals;
submaps_new.feats = initvals;
submaps_new.words = initvals;
submaps_new.ids = 1:numSubmaps;
% Rerank depth cannot be bigger than the number of submaps
if numSubmaps < opts.rerankDepth
    opts.rerankDepth = numSubmaps;
end
submaps_new.rerankDepth = opts.rerankDepth;

centers_used = [];
for i = 1:numSubmaps
    
    centerTravelled = true;
    center_ndx = 1;
    while centerTravelled
        center = DT.ConnectivityList(i,1);
        if ~ismember(center, centers_used)
            centerTravelled = false;
        else
            center_ndx = center_ndx + 1;
        end
    end
    sm_ndx = DT.ConnectivityList(i, center_ndx);
    sm_ctr = submaps.center{sm_ndx};
    
    submaps_new.center{i} = sm_ctr;
    submaps_new.H_to_world{i} = submaps.H_to_world{sm_ndx};
    submaps_new.feats{i} = submaps.feats{sm_ndx};
    submaps_new.words{i} = submaps.words{sm_ndx};
    
    for j = 1:2
        center_ndx = ndx(center_ndx+1, 3);
        sm_ndx2 = DT.ConnectivityList(i, center_ndx);
        sm_ctr2 = submaps.center{sm_ndx2};
        
        feats2 = submaps.feats{sm_ndx2};
        feats2(1, :) = feats2(1, :) + sm_ctr2(1) - sm_ctr(1);
        feats2(2, :) = feats2(2, :) + sm_ctr2(2) - sm_ctr(2);
        submaps_new.feats{i} = [submaps_new.feats{i}, feats2];
        
        words2 = submaps.words{sm_ndx2};
        submaps_new.words{i} = [submaps_new.words{i}, words2];
    end
    
    histogram = visualindex_get_histogram(model, submaps_new.words{i}(2,:));
    histograms(:, i) = full(histogram);
    
end

submaps_new.histograms = sparse(histograms);

if opts.showPlot
    plot_delaunay(DT); hold on
    plot(sm_ctrs(1,:),sm_ctrs(2,:),'kx',...
         'MarkerSize',12,'LineWidth',2)
    plot(sm_ctrs(1,:),sm_ctrs(2,:),'ko',...
         'MarkerSize',12,'LineWidth',2)
    set(gca, 'YDir', 'reverse')
    axis equal, axis tight
end

function index = ndx(index, numPoints)

if index > numPoints
    index = index - numPoints;
elseif index <= 0
    index = numPoints + index;
end