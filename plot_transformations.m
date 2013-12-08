function offsets = plot_transformations(cor, index, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 08/12/2013
% -------------------------------------------------------------------------
%
% PLOT_TRANSFORMATIONS
% plot_transformations(cor, index, varargin, 'LineColour', valLineColour,
% 'plotOnImage', valPlotOnImage)
%
% Plots the borders of all images in the global map. By using the 'hold on'
% command, this plot can be superimposed on other plots or images (e.g.
% those from plot_features or plot_feature_matches). This is done using the
% function plot_everything.
%
% Inputs:
%   - world:    World structure containing global features. Type 'help 
%               build_world' for more info
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
%
% Outputs:
%   - offset:   [xOffset; yOffset] values needed for superimposing the plot
%               on an image (see superimpose_plots for an example)

opts.LineColour = 'black';
opts.plotOnImage = false;
opts.dontPlot = false;
opts = vl_argparse(opts, varargin);

ims_mappable = find(cellfun(@(x)(~isempty(x)), cor.H_to_ref));
vertices = cell(1,length(ims_mappable));

xmin = 1; ymin = 1;
for i = 1:length(ims_mappable)
    img = ims_mappable(i);
    im_info = imfinfo(index.index.names{img});
    imsize = [im_info.Width; im_info.Height];
    [vertices{i}, xmin_tmp, ymin_tmp] = ...
        transform_rect(cor.H_to_ref{img}, imsize);
    % Update minimum values
    if xmin_tmp < xmin, xmin = xmin_tmp; end
    if ymin_tmp < ymin, ymin = ymin_tmp; end
end

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
    % pause(2); Add pause to visually check which image is which
end

% Reverse y-axis so that plot is aligned with image mosaic
set(gca, 'YDir', 'reverse')
axis equal, hold off

function [new_vertices, xmin, ymin] = ...
    transform_rect(H, imsize)

H_top = H(1:2,:);
H_bot = H(3,:);

% Vertices of image in homogeneous co-ordinates
im_vertices =   [  0    imsize(1)   imsize(1)           0   0
                   0            0   imsize(2)   imsize(2)   0
                   1            1           1           1   1   ];

% Transformed image
new_vertices = (H_top * im_vertices) ./ ...
    [H_bot * im_vertices; H_bot * im_vertices];

xmin = min(new_vertices(1,:));
ymin = min(new_vertices(2,:));