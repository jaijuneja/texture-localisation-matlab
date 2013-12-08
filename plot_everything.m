function plot_everything(world, cor, model, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 08/12/2013
% -------------------------------------------------------------------------
%
% PLOT_EVERYTHING
% plot_everything(world, cor, model, 'matchesOnly', valMatchesOnly,
% 'showMosaic', valShowMosaic, 'showImgBorders', valShowImgBorders)
%
% Superimposes plots of the image mosaic, features and image edge lines.
% Various optional input properties allow you to select which of these
% plots to include.
%
% Inputs:
%   - world:    World structure containing global features. Type 'help 
%               build_world' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%
%   Optional Properties:
%       - matchesOnly:      Set to true if you only want to show features
%                           that are matched between multiple images; false
%                           by default
%       - showMosaic:       Set to false if you don't want to display the
%                           image mosaic under the features; true by
%                           default
%       - showImgBorders:   Set to false if you don't want to display lines
%                           along the borders of images; true by default

opts.matchesOnly = false;
opts.showMosaic = true;
opts.showImgBorders = true;
opts = vl_argparse(opts, varargin);

% Get offset parameters
offsets = plot_transformations(cor, model, 'plotOnImage', true, ...
    'dontPlot', true);
xOffset = offsets(1); yOffset = offsets(2);

% Plot image mosaic
if opts.showMosaic
    mosaic = get_mosaic_pieces(model, cor);
    image_map = build_mosaic(model, mosaic, cor);
    imagesc(image_map); 
    hold on
end

% Plot features
if opts.matchesOnly
    plot_feature_matches(world, cor, 'xOffset', xOffset, 'yOffset', yOffset);
else
    plot_features(world, 'xOffset', xOffset, 'yOffset', yOffset);
end

hold on

% Plot image edges lines
if opts.showImgBorders
    plot_transformations(cor, model, 'plotOnImage', true, 'LineColour', 'g');
end

% Turn of axis if plotting mosaic
if opts.showMosaic, axis off, end
axis equal

end