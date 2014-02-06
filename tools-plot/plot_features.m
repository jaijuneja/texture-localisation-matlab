function plot_features(world, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% PLOT_FEATURES
% plot_features(world, 'plotStyle', varPlotStyle, 'numBins', varNumBins,
% 'numSamples', varNumSamples, 'xOffset', valXOffset, 'yOffset', valYOffset)
%
% Plots all of the features contained the the world structure in the global
% co-ordinate frame. Can either plot all individual points, a 2d histogram
% of points, or a random sample of points for very large feature maps. The
% histogram view is particularly appropriate for interpreting feature maps
% where large clusters of features exist.
%
% Inputs:
%   - world:    World structure containing global features. Type 'help 
%               build_world' for more info
%
% 	Optional Properties:
%       - plotStyle:    Can take values 'raw', 'histogram' or
%                       'randomSample'. Takes the value 'raw' by default
%       - numBins:      Only applicable for histogram plots. Has a default
%                       value of round(sqrt(world.num_features))
%       - numSamples:   Only applicable for sampled feature map plots
%                       Takes a default value of 1e5
%       - xOffset:      Offset in x-direction; needed for example when
%                       superimposing plot on an image
%       - yOffset:      Offset in y-direction
%       - scaleFactor:  Apply scale factor to feature co-ordinates
%       - showLegend:   Set to true by default

num_mappable_feats = sum(world.features_mappable);

opts.plotStyle = 'raw';
opts.numBins = round(sqrt(world.num_features));
opts.numSamples = 1e5;
opts.xOffset = 0;
opts.yOffset = 0;
opts.scaleFactor = 1;
opts.showLegend = true;
opts = vl_argparse(opts, varargin);

if ~isequal(opts.scaleFactor, 1)
    world = transform_world(world, opts.scaleFactor);
end

% If the number of samples specified is greater than the number of mappable
% features, there is no need to sample points. Just use raw plot method!
if strcmp(opts.plotStyle, 'randomSample') && opts.numSamples > num_mappable_feats
    opts.plotStyle = 'raw';
end

if strcmp(opts.plotStyle, 'raw')
    % Determine which features have been matched between multiple images
    matched = world.features_global(2,:) > 1;

    plot(world.features_global(3, ~matched & world.features_mappable) + opts.xOffset, ...
        world.features_global(4, ~matched & world.features_mappable) + opts.yOffset, 'b+', ... 
        world.features_global(3, matched) + opts.xOffset, ...
        world.features_global(4, matched) + opts.yOffset, 'r+');
    set(gca,'YDir','Reverse')
    axis equal
    if opts.showLegend, legend('Unmatched Features', 'Matched Features'); end
    title('Global Map of Features');
    
elseif strcmp(opts.plotStyle, 'histogram')

    % Collect the features that are mappable
    x_feats = world.features_global(3,world.features_mappable)';
    y_feats = world.features_global(4,world.features_mappable)';

    % Create evenly spaced bins in the x and y directions
    x_bins = linspace(min(x_feats), max(x_feats), opts.numBins);
    y_bins = linspace(min(y_feats), max(y_feats), opts.numBins);

    % Use nearest-neighbour interpolation to align feature positions with bins
    x_interp = interp1(x_bins, 1:numel(x_bins), x_feats, 'nearest')';
    y_interp = interp1(y_bins, 1:numel(y_bins), y_feats, 'nearest')';

    % Build a histogram of feature positions on the surface
    histo = accumarray([x_interp' y_interp'], 1, [opts.numBins opts.numBins]);

    surf(histo, 'LineStyle', 'none')

elseif strcmp(opts.plotStyle, 'randomSample')
    % Determine which features have been matched
    matched = world.features_global(2,:) > 1;

    % Remove un-mappable points
    matched = matched(world.features_mappable);
    
    % Randomly sample feature points and plot them
    sample_points = randsample(num_mappable_feats, opts.numSamples);
    feat_points = world.features_global(3:4, world.features_mappable);
    feat_points = feat_points(:, sample_points);
    matched = matched(sample_points);
    
    plot(feat_points(1, ~matched) + opts.xOffset, ...
        feat_points(2, ~matched) + opts.yOffset, 'b+', ... 
        feat_points(1, matched) + opts.xOffset, ...
        feat_points(2, matched) + opts.yOffset, 'r+');
    set(gca,'YDir','Reverse')
    axis equal
    legend('Unmatched Features', 'Matched Features');
    title('Global Map of Features');
else
    error('Property plotStyle can only take the values "raw", "histogram" and "randomSample"')
end

end