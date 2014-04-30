function points = set_true_pose(model, cor, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 01/04/2014
% -------------------------------------------------------------------------
%
% SET_TRUE_POSE
% [cor, H] = get_poses(model, cor, 'ScaleFactor', valScaleFactor)
%
% Computes the user-inputted 'true pose' of an image. The user selects four
% points that correspond to the corners of an image and the co-ordinates of
% these points in the world frame are calculated.
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%
% Outputs:
%   - points:   Points [x1 x2 x3 x4; y1 y2 y3 y4] on the world plane
%               corresponding to the corners of the image
% 
% See also: GET_POSES()

if isempty(cor.H_to_world)
    error('Global transformations first need to be set using get_poses')
end

opts.numPoses = 1;
opts.scaleFactor = 1;
opts = vl_argparse(opts, varargin);

if ~isequal(opts.scaleFactor, 1)
    cor = transform_world(cor, opts.scaleFactor);
end

cor.H_to_ref = cor.H_to_world;

% Plot the original mosaic
mosaic = get_mosaic_pieces(model, cor);
image_map = build_mosaic(model, mosaic, cor);
figure; imagesc(image_map);
hold on, axis equal, axis tight

points = cell(1, opts.numPoses)
for k = 1:opts.numPoses
	% Prompt user to input 4 points that should map to rectilinear co-ordinates
	numPoints = 4;
	x = zeros(1, numPoints);
	y = zeros(1, numPoints);

	for i = 1:numPoints
	    [x(i), y(i)] = ginput(1);
	    plot(x(i), y(i), 'r+'); hold on
	end

	plot([x x(1)], [y y(1)], 'r', 'LineWidth', 2); hold on
	drawnow;

	% Convert rectangle from mosaic co-ordinates to world co-ords (mosaic is in
	% positive co-ordinates only, whereas world includes negative values)
	offsets = plot_transformations(model, cor, 'dontPlot', true);
	x = x - offsets(1);
	y = y - offsets(2);
	pts = [x; y];

	points{k} = [x; y];
end

% Scale back to world if necessary
if ~isequal(opts.scaleFactor, 1)
    T_to_world = [1 0 0; 0 1 0; 0 0 opts.scaleFactor];
	for i = 1:opts.numPoses
    	points{i} = transform_points(points{i}, T_to_world);
    end
end

end