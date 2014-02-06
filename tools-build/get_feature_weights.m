function weights = get_feature_weights(world, cor, varargin)
opts.window_size = 50;

matched_feats = find(world.features_global(2,:) > 1);
feats_glob = world.features_global(3:4, matched_feats);

matching_scores = zeros(1, length(world.feature_map));
nearby_feats = matching_scores;

for i = matched_feats
    feat_glob = world.features_global(3:4, i);
    num_nearby_feats = feats_in_range(feat_glob, ...
        feats_glob, opts.window_size);

    feats_loc = world.feature_indices(:, i);
    feats_loc = feats_loc(feats_loc ~= 0);
    
    for j = 1:length(feats_loc)
        feat_loc = feats_loc(j);
        img_id = world.frames_local(2, feat_loc);
        
        matching_score = size(cell2mat(cor.feature_matches{img_id}), 2);
        matching_scores(feat_loc) = matching_score;

        nearby_feats(feat_loc) = num_nearby_feats;
    end
end

% Normalise
matching_scores = matching_scores / max(matching_scores);
weights = matching_scores ./ (nearby_feats + 1);
weights = weights / max(weights);

% ------------------------------------------------------
function num_feats = feats_in_range(feat, feats, window)
% ------------------------------------------------------
num_feats = ...
    (feats(1,:) > feat(1) - window) & ...
    (feats(1,:) < feat(1) + window) & ...
    (feats(2,:) > feat(2) - window) & ...
    (feats(2,:) < feat(2) + window);
num_feats = sum(num_feats);