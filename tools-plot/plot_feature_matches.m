function plot_feature_matches(world, cor, varargin)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 08/11/2013
% -------------------------------------------------------------------------
%
% PLOT_FEATURE_MATCHES
% plot_feature_matches(world, cor, 'xOffset', valXOffset, 'yOffset',
% valYOffset)
%
% Plots all global features that have multiple local matches across images.
% Also plots the local features that are matched to each global feature.
%
% Inputs:
%   - world:    World structure containing global features. Type 'help 
%               build_world' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%
%   Optional Properties:
%       - globalFeatsOnly:  Only show global features (hide local matches)
%                           when set to true; false by default
%       - xOffset:          Offset in x-direction; needed for example when
%                           superimposing plot on an image
%       - yOffset:          Offset in y-direction

opts.globalFeatsOnly = false;
opts.xOffset = 0;
opts.yOffset = 0;
opts = vl_argparse(opts, varargin);

% Find global features with multiple local matches
matched_global = world.features_global(2,:) > 1;
global_feats = world.features_global(3:4, matched_global);

% Plot global features in red
plot(global_feats(1,:) + opts.xOffset, global_feats(2,:) + opts.yOffset, 'r+');
hold on
labels{1} = 'Global Features';

local_feats = [];
if ~opts.globalFeatsOnly
    % Get indices of local features matched to the global ones
    matched_local = world.feature_indices(:, matched_global);
    matched_local = nonzeros(matched_local);
    local_frames = world.frames_local(:, matched_local); % Local SIFT ellipses

    % Transform local features to global frame
    ims_matched = find(cellfun(@(x)(~isempty(x)), cor.H_to_ref));
    local_feats = zeros(6, length(matched_local));
    for i = 1:length(ims_matched)
        img = ims_matched(i);
        img_feats = ismember(local_frames(2,:), img);
        local_feats(:, img_feats) = transform_frames(...
            local_frames(3:8, img_feats), cor.H_to_ref{img});
    end
        
    % Plot local features in blue. Plot should show the spread of local
    % features about each global feature in the global frame
    plot(local_feats(1,:) + opts.xOffset, local_feats(2,:) + opts.yOffset, 'bx');
    labels{2} = 'Local Features';
end

legend(labels);
title('Plot of Matched Features');
hold off

end