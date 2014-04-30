function cor = build_correspondence(model, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% BUILD_CORRESPONDENCE
% cor = get_correspondence(model, 'refImg', valRefImg, 'percentThresh',
% valPercentThresh, 'numThresh', valNumThresh)
%
% Generate a correspondence structure using the index of images.
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%
%   Optional Properties:
%       - 'refImg':         The ID of the reference image (which will
%                           coincide with the global co-ordinate frame)
%       - 'percentThresh':  The percentage of the maximum matching score
%                           above which a matched image must be to be
%                           considered a good match (this is a necessary
%                           but not sufficient condition). Set to 0.3 (30%)
%                           by default
%       - 'numThresh':      The value above which the matching score must
%                           be for an image to be considered a good match.
%                           Set to 20 by default
%
% Outputs:
%   - cor:	Correspondence structure containing links between different
%         	images (graph representation using an adjacency matrix):
%               *   cor.id is an ordered vector of image ids [1 2 ... n]
%               *   cor.ref_img is the ID of the image to be used as
%                   the initial reference co-ordinate frame
%               *   cor.img_matches is a 1xn cell array where each index
%                   cor.img_matches{i} contains the ids of the images
%                   matched to image i
%               *   cor.feature_matches is 1xn cell array where each index
%                   cor.feature_matches{i}{j} contains a 2xn vector of
%                   matched features between the jth matching image of
%                   query image i
%               *   cor.scores is a 1xn cell array where each index
%                   cor.scores{i}(j) is the matching score for the jth
%                   match to image i
%               *   cor.H is a 1xn cell array where each index cor.H{i}{j}
%                   is a 3x3 transformation matrix relating the jth matched
%                   image to image i
%               *   cor.H_to_ref is a 1xn cell array where each index
%                   cor.H_to_ref{i} is a 3x3 transformation matrix relating
%                   image i to the reference image (image 1 by default)
%               *   cor.H_world_toref is a 3x3 planar homography
%                   transforming world co-ordinates to pixel co-ordinates
%                   in the reference image up to an unknown scale. it is
%                   initialised as a blank array [], and is populated in
%                   the function get_poses
%               *   cor.H_to_world is similar to cor.H_to_ref, but relates
%                   each image to the world plane. In this function it is
%                   initialised as a blank array []. It is populated in
%                   the function get_poses
%               *   cor.adjacency is an nxn adjacency matrix where
%                   cor.adjacency(i,j) indicates a good match between
%                   images i and j. Ideally it should be symmetric but this
%                   may not be the case. A good way to view the matrix is
%                   by calling imagesc(cor.adjacency)
%               *   cor.graph is a biograph representation of the adjacency
%                   matrix. Use view(cor.graph) to view it
%               *   cor.intrinsics is the calibration matrix of the camera
%                   used. It is initialised as an empty array [] and is
%                   populated using the function set_intrinsics

% Score threshold is a percentage of the maximum score AND greater than a
% specified value (currently hardcoded and set to 20 below)
opts.percentThresh = 0.4;
opts.refImg = 1;
opts.numThresh = 20;
opts = vl_argparse(opts, varargin);

if opts.percentThresh > 1
    error('Property "percentThresh" must be less than 1')
end

num_imgs = length(model.index.ids);

% Initialise cor structure
fprintf('Initialising build... \n')
cor.id = zeros(1, num_imgs);
cor.img_matches = cell(1, num_imgs);
cor.feature_matches = cell(1, num_imgs);
cor.scores = cell(1, num_imgs);
cor.H = cell(1, num_imgs);
cor.H_to_ref = cell(1, num_imgs); % Transformation to reference image
cor.adjacency = zeros(num_imgs);

cor.H_world_toref = []; % Populate this using get_poses
cor.H_to_world = []; % Populate this using get_poses
cor.intrinsics = []; % Populate this using set_intrinsics

fprintf('Build progress: 000%%')
counter = 0;
for i = model.index.ids
    counter = counter + 1;
    [ids, scores, result] = visualindex_query(model, i, ...
        'geometricVerification', true);
    
    % Remove empty items
    nonempty = find(cellfun(@(x)(~isempty(x)), result.H));
    result.H = result.H(nonempty);
    result.ids = result.ids(nonempty);
    scores = scores(nonempty);
    result.matches = result.matches(nonempty);
    ids = ids(nonempty);
    
    % If there are no matches then skip to next iteration
    if length(scores) < 2, continue; end
    
    % Otherwise find the good matches
    score_thresh = scores(2) * opts.percentThresh; % of best match score    
    goodmatch = find(scores > score_thresh & scores > opts.numThresh ...
        & ids ~= result.query);
    cor.id(i) = result.query;
    cor.img_matches{i} = ids(goodmatch);
    cor.scores{i} = scores(goodmatch);
    cor.H{i} = result.H(goodmatch);
    cor.feature_matches{i} = result.matches(goodmatch);
    cor.adjacency(i, ids(goodmatch)) = 1;
    
    % Display percentage of build that is complete
    pc_complete = int64((counter / num_imgs) * 100);
    if pc_complete < 10, zer = '00'; 
    elseif pc_complete < 100, zer = '0';
    else zer = '';
    end
    fprintf(['\b\b\b\b' strcat(zer, num2str(pc_complete)) '%%'])
end

% Remove edges that are not bi-directional
[im_ids, match_ids] = find(cor.adjacency .* cor.adjacency' ~= cor.adjacency);
for i = 1:length(im_ids)
    match_ndx = (cor.img_matches{im_ids(i)} == match_ids(i));
    cor.img_matches{im_ids(i)}(match_ndx) = [];
    cor.feature_matches{im_ids(i)}(match_ndx) = [];
    cor.scores{im_ids(i)}(match_ndx) = [];
    cor.H{im_ids(i)}(match_ndx) = [];
    % Thus we zero out asymmetric terms
    cor.adjacency(im_ids(i), match_ids(i)) = 0;
end
fprintf('\n')

% Construct cell array of image names
node_names = cell(1, numel(model.index.ids));
node_names(:) = {'Image '};
node_ids = strtrim(cellstr(num2str(model.index.ids'))');
node_names(:) = strcat(node_names(:), node_ids(:));
% Construct biograph of images using adjacency matrix
cor.graph = biograph(tril(cor.adjacency), node_names, 'ShowArrows', 'off');

% Use view(graph) to display biograph
% Use imagesc(adjacency), axis equal to display adjacency matrix

% Set reference image and determine H_to_ref for all images
cor = set_refimg(cor, opts.refImg);
end