function world = visualindex_build_world(model, cor)

% Initialise global map. Add all features in first image to map
% First image acts as the reference frame
numWords = length(model.index.words{1});
world.features_global = 1:numWords;
world.features_local = cell(1, numWords); % Index of image and feature e.g. [1 3] = image 1, feature 3
world.words_global = cell(1, numWords);
world.frames_global = [];

for j = 1:numWords
    world.features_local{j} = [j; 1; j]; % [global_id; img_id; local_id]
    world.words_global{j} = [j; 1; model.index.words{1}(j)]; % [global_id; img_id; word]
end

% Incrementally add features from other images to global map
for i = 1:length(cor.id)
    
    numImgMatches = length(cor.img_matches{i});
    
    for j = 1:numImgMatches
        % NOTE: CURRENTLY ALLOWS FOR DOUBLE COUNTING OF FEATURES - e.g. IF
        % IMG 1 MATCHED TO IMG 2 THEN IMG 2 ALSO MATCHED TO IMG 1 LATER ON.
        % HENCE FEATURES CAN BE DOUBLE-COUNTED. DO WE ALLOW FOR THIS?
        matchedImgID = cor.img_matches{i}(j);
        im1_feats = cor.feature_matches{i}{j}(1,:);
        im2_feats = cor.feature_matches{i}{j}(2,:);
        
        initialised_features = cell2mat(world.features_local);
        % find features that can be matched to existing features in the map
        [matched_features, matched_ndx] = ismember(im1_feats, ...
            initialised_features(3,:) .* (initialised_features(2,:)==i));
        for k = 1:length(matched_features)
            im2_feat_toadd = im2_feats(k);
            % If the feature has been matched to one in global map...
            if matched_features(k) == 1
                
                global_feature_ndx = initialised_features(1,matched_ndx(k));
                
                % ... Then insert the feature to the correct location in
                % global map
                world.features_local{global_feature_ndx} = ...
                    [world.features_local{global_feature_ndx} ...
                    [global_feature_ndx; matchedImgID; im2_feat_toadd]];
                
                % Update frames and words in same way
                world.words_global{global_feature_ndx} = ...
                    [world.words_global{global_feature_ndx} ...
                    [global_feature_ndx; matchedImgID; ...
                    model.index.words{matchedImgID}(im2_feat_toadd)]];
                
            else
            % Otherwise, create a new feature in global map
            new_global_id = length(world.features_global)+1;
            world.features_global(new_global_id) = new_global_id;
            world.features_local{new_global_id} = ...
                [new_global_id; matchedImgID; im2_feat_toadd];
            
            % Update frames and words in same way
            world.words_global{new_global_id} = ...
                [new_global_id; matchedImgID; ...
                model.index.words{matchedImgID}(im2_feat_toadd)];
                
            end
        end
        
        % Add unmatched features to map
        num_features = size(model.index.frames{matchedImgID}, 2);
        feature_vector = 1:num_features;
        unmatched_features = ~ismember(feature_vector, im2_feats);
        unmatched_features = feature_vector(unmatched_features);
        numUnmatchedFeats = length(unmatched_features);
        
        % Add new global feature for each unmatched local feature
        new_global_id = length(world.features_global)+1;
        new_global_ids = new_global_id:(new_global_id + numUnmatchedFeats - 1);
        world.features_global = [world.features_global new_global_ids];
        
        % Add link to local feature for each unmatched local feature
        for k = 1:numUnmatchedFeats
            world.features_local{new_global_ids(k)} = ...
                [new_global_ids(k); matchedImgID; unmatched_features(k)];
            
            world.words_global{new_global_ids(k)} = ...
                [new_global_ids(k); matchedImgID; ...
                model.index.words{matchedImgID}(unmatched_features(k))];
        end
        
    end
    
end

% Need to work out the diff between frames and descriptors
% Might need to edit visualindex files to include exact descriptors in
% index
% At the moment using visual words instead of SIFT descriptors

end