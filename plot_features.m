% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% plot_features(world) plots all of the features contained the the world
% structure in the global co-ordinate frame.
%
% Inputs:
%   - world:    World structure containing global features. Type 'help 
%               build_world' for more info
%
function plot_features(world)
% Determine which features have been matched
histo = histc(world.feature_map(1,:), 1:world.num_features);
matched = histo > 1;

plot(world.frames_global(3, ~matched), world.frames_global(4, ~matched), 'b+', ... 
    world.frames_global(3, matched), world.frames_global(4, matched), 'r+', ...
    'LineStyle', 'none');
set(gca,'YDir','Reverse')
axis equal
legend('Unmatched Features', 'Matched Features');
title('Global Map of Features');
end