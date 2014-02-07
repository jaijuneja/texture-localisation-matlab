function cams = localize_image(img, model, cor, world, varargin)

opts.percentThresh = 0.5;
opts.numThresh = 20;
opts.intrinsicMat = cor.intrinsics;
opts.maxMatches = inf;
opts.plotPoses = true;
opts.showMosaic = false;
opts.scaleFactor = 1;
opts = vl_argparse(opts, varargin);

fprintf('Initialising localisation... \n');

if ~isequal(opts.scaleFactor, 1)
    cor = transform_world(cor, opts.scaleFactor);
end

% Set geometricVerification to false for coarse localisation approach
[ids, scores, result] = visualindex_query(model, img, ...
    'geometricVerification', true);
score_thresh = scores(1) * opts.percentThresh; % of best match score
goodmatch = find(scores > score_thresh & scores > opts.numThresh);
if length(goodmatch) > opts.maxMatches
    goodmatch = goodmatch(1:opts.maxMatches);
end
img_matches = ids(goodmatch);
scores = scores(goodmatch);
H = result.H(goodmatch);
feature_matches = result.matches(goodmatch);

% Find all possible poses
numPoses = length(goodmatch);
imsize = size(img);
imsize = [imsize(2); imsize(1)];
cams.imsize = imsize;

for i = 1:numPoses
    match = img_matches(i);
    H_match = cor.H_to_world{match};
    if isempty(H_match)
        continue;
    end
    H_to_match = H{i};
    H_to_world = H_match / H_to_match;
    H_world_toimg = inv(H_to_world);
    H_world_toimg = H_world_toimg / H_world_toimg(3,3);
    [cams.R{i}, cams.t{i}] = ...
        decompose_homog(H_world_toimg, opts.intrinsicMat);
    cams.H_to_world{i} = H_to_world;
    cams.feature_matches{i} = feature_matches{i};
    cams.scores(i) = scores(i);
end

numPoses = length(cams.R);
pose_str = 'poses';
if isequal(numPoses, 1), pose_str = 'pose'; end
fprintf(['Image localised, ' num2str(numPoses) ' ' pose_str ' found. \n'])

if opts.plotPoses
    plot3d_poses(model, cor, world, 'queryCameras', cams, ...
        'onlyDrawQueryCams', true, 'showMosaic', opts.showMosaic)
end

% Now we have coarse localisation of image. Use this to determine which
% submap(s) the camera falls into

% Using feature matching with submaps, get finer localisation
end