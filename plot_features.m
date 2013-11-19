function plot_features(world)
% Determine which features have been matched
histo = histc(world.feature_map(1,:), 1:world.num_features);
matched = histo > 1;

figure;
plot(world.frames_global(3, ~matched), world.frames_global(4, ~matched), 'b+', ... 
    world.frames_global(3, matched), world.frames_global(4, matched), 'r+', ...
    'LineStyle', 'none');
set(gca,'YDir','Reverse')
axis equal
legend('Unmatched Features', 'Matched Features');
title('Global Map of Features');
end