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
%                   the reference co-ordinate frame
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
%               *   cor.adjacency is an nxn adjacency matrix where
%                   cor.adjacency(i,j) indicates a good match between
%                   images i and j. Ideally it should be symmetric but this
%                   may not be the case. A good way to view the matrix is
%                   by calling imagesc(cor.adjacency)
%               *   cor.graph is a biograph representation of the adjacency
%                   matrix. Use view(cor.graph) to view it

% Score threshold is a percentage of the maximum score AND greater than a
% specified value (currently hardcoded and set to 20 below)
opts.percentThresh = 0.5;
opts.refImg = 1;
opts.numThresh = 20;
opts = vl_argparse(opts, varargin);

if opts.percentThresh > 1
    error('Property "percentThresh" must be less than 1')
end

% Initialise cor structure
cor.id = zeros(1, length(model.index.ids));
cor.img_matches = cell(1, length(model.index.ids));
cor.feature_matches = cell(1, length(model.index.ids));
cor.scores = cell(1, length(model.index.ids));
cor.H = cell(1, length(model.index.ids));
cor.H_to_ref = cell(1, length(model.index.ids)); % Transformation to reference image
cor.adjacency = zeros(length(model.index.ids));

for i = model.index.ids
    [ids, scores, result] = visualindex_query(model, i);
    score_thresh = scores(2) * opts.percentThresh; % 60% of best score by default
    goodmatch = find(scores > score_thresh & scores > opts.numThresh ...
        & ids ~= result.query);
    cor.id(i) = result.query;
    cor.img_matches{i} = ids(goodmatch);
    cor.scores{i} = scores(goodmatch);
    cor.H{i} = result.H(goodmatch);
    cor.feature_matches{i} = result.matches(goodmatch);
    cor.adjacency(i, ids(goodmatch)) = 1;
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
    
% Construct cell array of image names
node_names = cell(1, numel(model.index.ids));
node_names(:) = {'Image '};
node_ids = strtrim(cellstr(num2str(model.index.ids'))');
node_names(:) = strcat(node_names(:), node_ids(:));
% Construct biograph of images using adjacency matrix
cor.graph = biograph(cor.adjacency, node_names);
% Use view(graph) to display biograph
% Use imagesc(adjacency), axis equal to display adjacency matrix

% Set reference image and determine H_to_ref for all images
cor = set_refimg(cor, opts.refImg);

end