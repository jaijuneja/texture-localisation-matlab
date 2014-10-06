function [matches, cams] = localisation_submaps(model, world, ...
    cor, img, varargin)

opts.submapStyle = 'kmeans';
opts.percentThresh = 0.5;
opts.numThresh = 15;
opts.maxMatches = inf;
opts.intrinsicMat = cor.intrinsics;
opts.plotPoses = true;
opts.showMosaic = false;
opts.scaleFactor = 1;
opts.numCenters = 20;
opts.featsPerSubmap = [];
opts.autoGenerateSubmaps = false;
opts.showSubmaps = false;
opts.truePosition = [];
opts.timeLocalisation = false;
opts = vl_argparse(opts, varargin);

[~, ref_dims] = plot_transformations(model, cor, 'dontPlot', 'true', ...
    'fromFrame', 'ref');
[~, wld_dims] = plot_transformations(model, cor, 'dontPlot', 'true');
scalediff = norm(ref_dims)/norm(wld_dims);

geomVerTol1 = 35 / scalediff;
geomVerTol2 = 30 / scalediff;
geomVerTol3 = 25 / scalediff;

% extract the features, visual words, and histogram for the query images
fprintf('Extracting and quantising features from query image... ')
if numel(img) > 1
    % im is an image; extract features
    [frames, descrs] = visualindex_get_features(model, img) ;
    words = visualindex_get_words(model, descrs) ;
else
    % im is an index in the database
    frames = model.index.frames{img} ;
    words = model.index.words{img} ;
end

fprintf('Done. \nBuilding submaps... ')
if opts.autoGenerateSubmaps
    opts.featsPerSubmap = size(words, 2) * 1.5;
end

if strcmp(opts.submapStyle, 'kmeans')
    submaps = build_submaps_kmeans(model, world, 'showPlot', opts.showSubmaps, ...
        'numCenters', opts.numCenters, 'featsPerSubmap', opts.featsPerSubmap);
elseif strcmp(opts.submapStyle, 'rect')
    submaps = build_submaps_rect(model, world, 'showPlot', opts.showSubmaps, ...
        'numSubmaps', opts.numCenters, 'featsPerSubmap', opts.featsPerSubmap);
else
    error('Property "submapStyle" can only take the values "kmeans" or "rect"')
end
fprintf(['Done. ' num2str(length(submaps.center)) ...
    ' submap(s) generated. \nInitialising localisation... '])

if opts.timeLocalisation
    tic
end

% First compute matching scores against submaps using histograms of words
histogram = visualindex_get_histogram(model, words);

% Compute histogram-based score
scores = histogram' * submaps.histograms;
[scores, perm] = sort(scores, 'descend');

depth = submaps.rerankDepth;
scores = scores(1:depth);
perm = perm(1:depth);
ids = submaps.ids(perm);

matches = cell(1,depth);
H = cell(1,depth);

fprintf('Done. \nLocalising image using geometric verification... \n')
% Update top scores using geometric verification
words2 = submaps.words(perm);
frames2 = submaps.feats(perm);

ambig_matches = zeros(1, depth);

% Initialise parallel processing
if matlabpool('size') == 0 % Check if pool is already open
    matlabpool open 2
end

parfor t = 1:depth
    % Find the features that are mapped to the same visual words

    % Iterate through all feature matches including those where ambiguities
    % exist (i.e. mapping to multiple visual words)
    [feats2, feats1] = ismember(words2{t}(2,:), words);
    m1 = feats1(feats1~=0);
    m2 = find(feats2);
    matches{t} = [m1(:)'; m2(:)'];
    ambig_matches(t) = length(m2)/length(unique(m1));
    
    % Alternatively, only iterate through first match where ambiguities exist
    % [~,m1,m2] = intersect(words,words2{t}(2,:));
    % matches{t} = [m1(:)';m2(:)'];
    
    [inliers, H{t}] = geometricVerification(frames, frames2{t}, matches{t}, ...
        'tolerance1', geomVerTol1, 'tolerance2', geomVerTol2, ...
        'tolerance3', geomVerTol3);
    
    if numel(inliers) >= 6
        scores(t) = scores(t) + numel(inliers);
    end
    matches{t} = matches{t}(:, inliers);
end

ambig_matches = mean(ambig_matches);
fprintf(['On average ' num2str(ambig_matches) ...
    ' matches per feature during geometric verification. \n']);

% rerank
[scores, perm] = sort(scores, 'descend') ;

score_thresh = scores(1) * opts.percentThresh; % of best match score
goodmatch = find(scores > score_thresh & scores > opts.numThresh);
if length(goodmatch) > opts.maxMatches
    goodmatch = goodmatch(1:opts.maxMatches);
end
scores = scores(goodmatch);
perm = perm(goodmatch);

ids = ids(perm) ;
matches = matches(perm) ;
H = H(perm) ;

% Find all possible poses
numPoses = length(goodmatch);
imsize = size(img);
imsize = [imsize(2); imsize(1)];
cams.imsize = imsize;

for i = 1:numPoses
    match = ids(i);
    H_match = submaps.H_to_world{match};
    
    H_to_match = H{i};
    H_to_world =  H_match / H_to_match;
    
    H_world_toimg = inv(H_to_world);
    H_world_toimg = H_world_toimg / H_world_toimg(3,3);

    [cams.R{i}, cams.t{i}] = ...
        decompose_homog(H_world_toimg, opts.intrinsicMat);
    cams.H_to_world{i} = H_to_world;
    cams.feature_matches{i} = matches{i};
    cams.scores(i) = scores(i);
end

numPoses = length(cams.R);
pose_str = 'poses';
if isequal(numPoses, 1), pose_str = 'pose'; end
fprintf(['Image localised, ' num2str(numPoses) ' ' pose_str ' found. \n'])

if opts.timeLocalisation
    toc
end

if opts.plotPoses
    if ~isempty(opts.truePosition), overlap = true; else overlap = false; end
    plot3d_poses(model, cor, world, 'queryCameras', cams, ...
        'onlyDrawQueryCams', true, 'showMosaic', opts.showMosaic, ...
        'truePosition', opts.truePosition, 'computeOverlap', overlap)
end

end