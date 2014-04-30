function submaps = build_submaps_rect(model, world, varargin)

opts.numCentersX = 10;
opts.numCentersY = 10;
opts.numSubmaps = [];
opts.featsPerSubmap = [];
opts.rerankDepth = 5;
opts.showPlot = true;
opts = vl_argparse(opts, varargin);

feats_glob = world.features_global(3:end, :);
num_feats = length(feats_glob);

if ~isempty(opts.featsPerSubmap)
    opts.numSubmaps = ceil(num_feats / opts.featsPerSubmap);
end

if ~isempty(opts.numSubmaps)
    xy_centers = ceil(sqrt(opts.numSubmaps));
    opts.numCentersX = xy_centers + 1;
    opts.numCentersY = xy_centers + 1;
end

% Get submap centres and indexes of features belonging to each submap
[sm_hist, sm_ctrs, sm_idxs] = histc_2d(feats_glob(1:2,:), ...
    opts.numCentersX, opts.numCentersY);
num_submaps = size(sm_ctrs, 2);

% Create matrix of colours
ColOrd = get(gca,'ColorOrder');
[m,~] = size(ColOrd);

for i = 1:num_submaps
    ColRow = rem(i,m);
    if ColRow == 0
      ColRow = m;
    end
    % Get the color
    Col = ColOrd(ColRow,:);

    if opts.showPlot
        plot(feats_glob(1, sm_idxs==i), feats_glob(2, sm_idxs==i), '.', ...
            'Color', Col, 'MarkerSize',12)
        hold on
    end
end

% Merge submaps so that they overlap
[sm_ctrs, sm_idxs] = merge_submaps(sm_hist, sm_ctrs, sm_idxs);
num_submaps = size(sm_ctrs, 2);

% Plot the submap centers
if opts.showPlot
    plot(sm_ctrs(1,:),sm_ctrs(2,:),'kx',...
         'MarkerSize',20,'LineWidth',2); hold on
    plot(sm_ctrs(1,:),sm_ctrs(2,:),'ko',...
         'MarkerSize',20,'LineWidth',2)
    set(gca, 'YDir', 'reverse')
    axis equal, axis tight
end

initvals = cell(1, num_submaps);

submaps.center = initvals;
submaps.H_to_world = initvals;
submaps.feats = initvals;
submaps.words = initvals;
submaps.ids = 1:num_submaps;

num_words = size(model.index.histograms, 1);
histograms = zeros(num_words, num_submaps);
% Rerank depth cannot be bigger than the number of submaps
if num_submaps < opts.rerankDepth
    opts.rerankDepth = num_submaps;
end
submaps.rerankDepth = opts.rerankDepth;

for i = 1:num_submaps
    submaps.center{i} = sm_ctrs(:, i);
    submaps.H_to_world{i} = [1 0 sm_ctrs(1, i); 0 1 sm_ctrs(2, i); 0 0 1];
    feats = world.features_global(3:end, sm_idxs{i});
    feats(1, :) = feats(1, :) - sm_ctrs(1, i);
    feats(2, :) = feats(2, :) - sm_ctrs(2, i);
    submaps.feats{i} = feats;
    submaps.words{i} = world.words_global(:, sm_idxs{i});
    
    histogram = visualindex_get_histogram(model, submaps.words{i}(2,:));
    histograms(:, i) = full(histogram);  
end

submaps.histograms = sparse(histograms);

function [new_ctrs, new_bins] = merge_submaps(hist, ctrs, bins)

[y_ctrs, x_ctrs] = size(hist);

num_new_ctrs = (x_ctrs-1)*(y_ctrs-1);
new_ctrs = zeros(2, num_new_ctrs);
new_bins = cell(1, num_new_ctrs);
for i = 1:x_ctrs-1
    for j = 1:y_ctrs - 1
        ndxs = [convert_ndx(i, j, x_ctrs), convert_ndx(i+1, j, x_ctrs), ...
            convert_ndx(i, j+1, x_ctrs), convert_ndx(i+1, j+1, x_ctrs)];
        new_ndx = convert_ndx(i, j, x_ctrs-1);
        new_ctrs(:, new_ndx) = [mean(ctrs(1, ndxs)); mean(ctrs(2, ndxs))];
        bin_ndxs = ismember(bins, ndxs);
        new_bins{new_ndx} = find(bin_ndxs);
    end
end

function [n, ctrs, bins] = histc_2d(pts, num_xctrs, num_yctrs)

xmin = min(pts(1,:));
xmax = max(pts(1,:));
ymin = min(pts(2,:));
ymax = max(pts(2,:));
dx_ctrs = (xmax-xmin)/num_xctrs;
x_ctrs = linspace(xmin + dx_ctrs/2, xmax - dx_ctrs/2, num_xctrs);
dy_ctrs = (ymax-ymin)/num_yctrs;
y_ctrs = linspace(ymin + dy_ctrs/2, ymax - dy_ctrs/2, num_yctrs);

x_vals = pts(1,:);
y_vals = pts(2,:);

dx_ctrs = diff(x_ctrs); % Bin width vector
dy_ctrs = diff(y_ctrs); % Bin width vector
assert(isequal(num_xctrs, size(x_ctrs,2)));
assert(isequal(num_yctrs, size(y_ctrs,2)));
edge_x = [-inf , x_ctrs(1:end-1)+dx_ctrs/2 , inf]; % Shifted edges for histc
edge_y = [-inf , y_ctrs(1:end-1)+dy_ctrs/2 , inf]; % Shifted edges for histc

ctrs = zeros(2, num_xctrs * num_yctrs);
ctrs(1,:) = repmat(x_ctrs, 1, num_yctrs);
for i = 1:num_yctrs
    ctrs(2, num_xctrs * (i - 1) + 1 : num_xctrs * i) = ...
        repmat(y_ctrs(i), 1, num_xctrs);
end

num_pts = size(pts, 2);
bins = zeros(1, num_pts);
n = zeros(num_xctrs+1,num_yctrs+1);
[n_x,bin_x] = histc(x_vals,edge_x); % Sort y1 into bins
for i = 1:num_xctrs
    if ~isequal(n_x, 0)
        ndxs = (bin_x == i);
        [n(i,:), bin_yi] = histc(y_vals(ndxs)',edge_y);
        bin_yi = convert_ndxs(i, bin_yi, num_xctrs);
        bins(ndxs) = bin_yi;
    end
end
n(:, end) = []; n(end, :) = [];

function y_ndxs = convert_ndxs(x_ndx, y_ndxs, xdim)

for i = 1:length(y_ndxs)
    y_ndxs(i) = convert_ndx(x_ndx, y_ndxs(i), xdim);
end

function ndx = convert_ndx(x_ndx, y_ndx, xdim)

ndx = xdim * (y_ndx - 1) + x_ndx;