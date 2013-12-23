function world = build_world(model, cor)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
%
% BUILD_WORLD
% world = build_world(model, cor)
% 
% Build world model comprised of global features linked to local features
% in different images.
%
% Inputs:
%   - model:    Index of images from visualindex. Type 'help
%               visualindex_build' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation using an adjacency matrix).
%               Type 'help build_correspondence' for more info
%
% Outputs:
%   - world:    World map structure:
%                   *   world.feature_map is a 3xn matrix for n total local
%                       features, each of which is mapped to a global
%                       feature in the form:
%                       [global_feat_id; img_id; local_feat_id]
%                   *   world.features_global is a 3xm matrix for m unique
%                       global features, containing estimates of each
%                       feature's global position. It takes the form:
%                       [global_feat_id; num_matched_feats; x_pos; y_pos]
%                   *   world.feature_indices is a sparse matrix with m
%                       columns (for m unique global features). Non-zero
%                       elements in the mth column are the indices of local
%                       features in world.feature_map which map to global
%                       feature m
%                   *   world.num_features = length(world.features_global)
%                   *   world.words_global is a 3xn matrix with the form:
%                       [global_feat_id; img_id; visual_word]
%                   *   world.frames_local is a 8xn matrix with the form:
%                       [global_feat_id; img_id; x_pos; y_pos; ...
%                        R(1,1); R(2,1); R(1,2); R(2,2)]
%                       where R a 2x2 matrix that transforms the unit
%                       circle to an orientated ellipse. Frames are defined
%                       in the local co-ordinates of img_id
%                   *   world.features_mappable is a 1xm logical array that
%                       is true for features whose global position is knwon
%                       and false otherwise

% Initialise global map. Add all features in first image to map.
% First image acts as the reference frame so we can populate feature info
% from the first image
world.num_features = length(model.index.words{cor.ref_img});

world.feature_map = [ 1:world.num_features
                      repmat(cor.ref_img, 1,world.num_features)
                      1:world.num_features    ];

frames_glob = transform_frames(model.index.frames{cor.ref_img}, ...
    cor.H_to_ref{cor.ref_img});

world.frames_local = [ 1:world.num_features
                        repmat(cor.ref_img, 1,world.num_features)
                        model.index.frames{cor.ref_img}   ];
                    
world.features_global = [ 1:world.num_features
                          ones(1,world.num_features)
                          frames_glob(1:2,:) ];
                      
world.feature_indices = 1:world.num_features;

world.feature_indices = sparse(world.feature_indices);

% At the moment using visual words instead of SIFT descriptors
world.words_global = [ 1:world.num_features
                       repmat(cor.ref_img, 1,world.num_features)
                       model.index.words{cor.ref_img}    ];
                    
world.features_mappable = true(1, world.num_features);

% Keep a tab of the number and IDs of images mapped
images_mapped = zeros(1, length(cor.img_order));
images_mapped(1) = cor.ref_img;
num_images_mapped = 1;

num_img_matches = length(cell2mat(cor.img_matches));
match_ids = zeros(1, num_img_matches);
num_ims_matched = 0;

% Incrementally add features from other images to global map
for i = cor.img_order
    
    numImgMatches = length(cor.img_matches{i});
    
    for j = 1:numImgMatches
        
        matchedImgID = cor.img_matches{i}(j);
        num_ims_matched = num_ims_matched + 1;

        % Checked whether the current image match has already been mapped
        % First create match_id by concatenating image IDs
        if i < matchedImgID
            im1_id = i;
            im2_id = matchedImgID;
        else
            im1_id = matchedImgID;
            im2_id = i;
        end
        
        match_id = str2double(strcat(num2str(im1_id),num2str(im2_id)));

        % If the current image match has already been mapped then skip
        % to next iteration
        if ismember(match_id, match_ids)
            match_ids(num_ims_matched) = match_id;
            continue;
        end
        
        match_ids(num_ims_matched) = match_id;
        
        im1_feats = cor.feature_matches{i}{j}(1,:);
        im2_feats = cor.feature_matches{i}{j}(2,:);
        
        % Find features that can be matched to existing features in the map
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
                new_global_id = world.num_features + 1;

                world = update_features(world, cor, model, ...
                    new_global_id, matchedImgID, im2_feat_toadd);
            end
        end
        
        % If the unmatched features from the matched image have already
        % been mapped, then skip to the next iteration. This step ensures
        % that unmatched features are not double counted
        
        if ismember(matchedImgID, images_mapped)
            continue;
        end
        
        num_images_mapped = num_images_mapped + 1;
        images_mapped(num_images_mapped) = matchedImgID;
        
        % Add unmatched features to map
        num_features = size(model.index.frames{matchedImgID}, 2);
        feature_vector = 1:num_features;
        unmatched_features = ~ismember(feature_vector, im2_feats);
        unmatched_features = feature_vector(unmatched_features);
        numUnmatchedFeats = length(unmatched_features);
        
        % Add a new global feature for each unmatched local feature
        new_global_id = world.num_features + 1;
        new_global_ids = new_global_id : world.num_features + numUnmatchedFeats;
        
        % Add link to local feature for each unmatched local feature
        for k = 1:numUnmatchedFeats
            world = update_features(world, cor, model, ...
                new_global_ids(k), matchedImgID, unmatched_features(k));             
        end
        
    end
    
end

end