function [offsets, dims] = plot_transformations(model, cor, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 08/12/2013
% -------------------------------------------------------------------------
%
% PLOT_TRANSFORMATIONS
% plot_transformations(model, cor, varargin, 'LineColour', valLineColour,
% 'plotOnImage', valPlotOnImage)
%
% Plots the borders of all images in the global map. By using the 'hold on'
% command, this plot can be superimposed on other plots or images (e.g.
% those from plot_features or plot_feature_matches). This is done using the
% function plot_everything.
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%
%   Optional Properties:
%       - LineColour:   Colour of the border lines
%       - plotOnImage:  Set to true when the plot needs to be superimposed
%                       on an image mosaic. Set to false by default
%       - dontPlot:     Set to true when you don't want to plot the result,
%                       but just need the offset values
%       - addDelay:     When set to true, adds a 1 second delay between
%                       plotting each image
%       - fromFrame:    Which reference frame to plot the world from.
%                       Either the world frame 'w', or the ref frame 'ref'
%
% Outputs:
%   - offset:   [xOffset; yOffset] values needed for superimposing the plot
%               on an image (see plot_everything for an example)
%   - dims:     [width; height] dimensions of the world

opts.LineColour = 'black';
opts.plotOnImage = false;
opts.dontPlot = false;
opts.addDelay = false;
opts.fromFrame = 'w';
opts = vl_argparse(opts, varargin);

ims_mappable = find(cellfun(@(x)(~isempty(x)), cor.H_to_ref));
vertices = cell(1,length(ims_mappable));

for i = 1:length(ims_mappable)
    img = ims_mappable(i);
    im_info = imfinfo(model.index.names{img});
    imsize = [im_info.Width; im_info.Height];
    if strcmp(opts.fromFrame, 'w')
        vertices{i} = transform_rect(cor.H_to_world{img}, imsize);
    else
        vertices{i} = transform_rect(cor.H_to_ref{img}, imsize);
    end
end

% Get offset values
vertices_mat = cell2mat(vertices);
xmin = min(vertices_mat(1,:));
ymin = min(vertices_mat(2,:));
xmax = max(vertices_mat(1,:));
ymax = max(vertices_mat(2,:));

dims = [xmax - xmin; ymax - ymin];

if xmin > 0, xmin = 0; end
if ymin > 0, ymin = 0; end
xOffset = abs(xmin) + 0.5;
yOffset = abs(ymin) + 0.5;

% Output offset values
offsets = [xOffset; yOffset];
% If we don't want to plot anything, skip the rest of the function
if opts.dontPlot
    return;
end

for i = 1:length(ims_mappable)
    % If plotting on an image, then adjust plot so that the minimum
    % x and y co-ordinate is 0.5
    if opts.plotOnImage
        vertices{i}(1,:) = vertices{i}(1,:) + xOffset;
        vertices{i}(2,:) = vertices{i}(2,:) + yOffset;
    end
    plot(vertices{i}(1,:), vertices{i}(2,:), opts.LineColour, ...
        'LineWidth', 1.5);
    hold on
    if opts.addDelay
        pause(1); % Add pause to visually check which image is which
    end
end

% Reverse y-axis so that plot is aligned with image mosaic
set(gca, 'YDir', 'reverse')
axis equal, hold off


% -----------------------------------------------
function new_vertices = transform_rect(H, imsize)
% -----------------------------------------------
% Transforms the corners of a rectangle with dimensions imsize(2) x
% imsize(1) by the transformation matrix H. Returns the new vertices.

H_top = H(1:2,:);
H_bot = H(3,:);

% Vertices of image in homogeneous co-ordinates (note that 5th vertex is
% same as first, so that when plotting them they form a closed polygon)
im_vertices =   [  1    imsize(1)   imsize(1)           1   1
                   1            1   imsize(2)   imsize(2)   1
                   1            1           1           1   1   ];

% Vertices of transformed image
new_vertices = (H_top * im_vertices) ./ ...
    [H_bot * im_vertices; H_bot * im_vertices];