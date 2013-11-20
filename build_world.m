% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% Build world model comprised of global features linked to local features
% in different images. Call world = build_world(model, cor).
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correlation structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help get_correlation' for more info
%
% Outputs:
%   - world:    World map structure:
%                   *   world.feature_map is a 3xn matrix for n global
%                       features, with the form:
%                       [global_feat_id; img_id; local_feat_id]
%                   *   world.num_features = length(world.feature_map)
%                   *   world.words_global is a 3xn matrix for n global
%                       features, with the form:
%                       [global_feat_id; img_id; visual_word]
%                   *   world.frames_global is a 8xn matrix for n global
%                       features, with the form:
%                       [global_feat_id; img_id; x_pos; y_pos; ...
%                        R(1,1); R(2,1); R(1,2); R(2,2)]
%                       where R a 2x2 matrix that transforms the unit
%                       circle to an orientated ellipse
%
function world = build_world(model, cor)

% Initialise global map. Add all features in first image to map
% First image acts as the reference frame
world.num_features  = length(model.index.words{1});
world.feature_map   = zeros(3, world.num_features);
world.words_global  = zeros(3, world.num_features);
world.frames_global = zeros(8, world.num_features);

world.feature_map   = [ 1:world.num_features
                        ones(1,world.num_features)
                        1:world.num_features    ];

% At the moment using visual words instead of SIFT descriptors
world.words_global  = [ 1:world.num_features
                        ones(1,world.num_features)
                        model.index.words{1}    ];
                    
world.frames_global = [ 1:world.num_features
                        ones(1,world.num_features)
                        model.index.frames{1}   ];

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
        
        % find features that can be matched to existing features in the map
        [matched_features, matched_ndx] = ismember(im1_feats, ...
            world.feature_map(3,:) .* (world.feature_map(2,:)==i));
        
        for k = 1:length(matched_features)
            im2_feat_toadd = im2_feats(k);
            
            % If the feature has been matched to one in global map
            if matched_features(k) == 1
                
                global_feature_id = world.feature_map(1,matched_ndx(k));
                
                world = update_features(world, cor, model, ...
                    global_feature_id, matchedImgID, im2_feat_toadd);
                
            else
                % Otherwise, create a new feature in global map
                world.num_features = world.num_features + 1;
                new_global_id = world.num_features;

                world = update_features(world, cor, model, ...
                    new_global_id, matchedImgID, im2_feat_toadd);                
            end
        end
        
        % Add unmatched features to map
        num_features = size(model.index.frames{matchedImgID}, 2);
        feature_vector = 1:num_features;
        unmatched_features = ~ismember(feature_vector, im2_feats);
        unmatched_features = feature_vector(unmatched_features);
        numUnmatchedFeats = length(unmatched_features);
        
        % Add a new global feature for each unmatched local feature
        new_global_id = world.num_features + 1;
        world.num_features = world.num_features + numUnmatchedFeats;
        new_global_ids = new_global_id : world.num_features;
        
        % Add link to local feature for each unmatched local feature
        for k = 1:numUnmatchedFeats
            world = update_features(world, cor, model, ...
                new_global_ids(k), matchedImgID, unmatched_features(k));             
        end
        
    end
    
end

end