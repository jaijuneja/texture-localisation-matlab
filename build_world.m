function world = build_world(model, cor)

% Initialise global map. Add all features in first image to map
% First image acts as the reference frame
world.num_features  = length(model.index.words{1});
world.feature_map   = zeros(3, world.num_features); % Index of image and feature e.g. [1 3] = image 1, feature 3
world.words_global  = zeros(3, world.num_features);
world.frames_global = zeros(8, world.num_features);

world.feature_map   = [ 1:world.num_features
                        ones(1,world.num_features)
                        1:world.num_features    ];
                    
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
        m = matchedImgID; % Shorthand used in some parts of code
        
        im1_feats = cor.feature_matches{i}{j}(1,:);
        im2_feats = cor.feature_matches{i}{j}(2,:);
        
        % find features that can be matched to existing features in the map
        [matched_features, matched_ndx] = ismember(im1_feats, ...
            world.feature_map(3,:) .* (world.feature_map(2,:)==i));
        
        for k = 1:length(matched_features)
            im2_feat_toadd = im2_feats(k);
            fndx = im2_feat_toadd; % Shorthand used in some parts of code
            
            % If the feature has been matched to one in global map...
            if matched_features(k) == 1
                
                global_feature_id = world.feature_map(1,matched_ndx(k));
                
                % ... Then insert the feature to the correct location in
                % global map
                world.feature_map(:, end+1) = ...
                    [global_feature_id; matchedImgID; im2_feat_toadd];
                
                % Update words in same way
                world.words_global(:, end+1) = ...
                    [global_feature_id; matchedImgID; ...
                    model.index.words{matchedImgID}(im2_feat_toadd)];
                
                if ~isempty(cor.H_to_ref{matchedImgID})
                    % Transform local frames to global co-ordinate system
                    % First obtain matrix transforming unit circle to
                    % orientated elliptical frame in local co-ordinates
                    T = [ ...
                        [model.index.frames{m}(3,fndx); model.index.frames{m}(4,fndx); 0] ...
                        [model.index.frames{m}(5,fndx); model.index.frames{m}(6,fndx); 0] ...
                        [model.index.frames{m}(1:2,fndx); 1] ...
                        ];
                    % Then pre-multiply this by the transformation to global
                    % co-ordinates to get the global feature ellipse
                    global_frame = cor.H_to_ref{matchedImgID} \ T;
                    global_frame = global_frame / global_frame(3,3);

                    world.frames_global(:, end+1) = ...
                        [
                        global_feature_id;
                        matchedImgID;
                        global_frame(1:2,3)
                        global_frame(1,1)
                        global_frame(2,1)
                        global_frame(1,2)
                        global_frame(2,2)                    
                        ];
                else
                    world.frames_global(:, end+1) = ...
                        [global_feature_id; matchedImgID; zeros(6,1)];
                end
                
            else
                % Otherwise, create a new feature in global map
                world.num_features = world.num_features + 1;
                new_global_id = world.num_features;

                world.feature_map(:, end+1) = ...
                    [new_global_id; matchedImgID; im2_feat_toadd];

                % Update words in same way
                world.words_global(:, end+1) = ...
                    [new_global_id; matchedImgID; ...
                    model.index.words{matchedImgID}(im2_feat_toadd)];

                if ~isempty(cor.H_to_ref{matchedImgID})
                    % Transform local frames to global co-ordinate system
                    % First obtain matrix transforming unit circle to
                    % orientated elliptical frame in local co-ordinates
                    T = [ ...
                        [model.index.frames{m}(3,fndx); model.index.frames{m}(4,fndx); 0] ...
                        [model.index.frames{m}(5,fndx); model.index.frames{m}(6,fndx); 0] ...
                        [model.index.frames{m}(1:2,fndx); 1] ...
                        ];
                    % Then pre-multiply this by the transformation to global
                    % co-ordinates to get the global feature ellipse
                    global_frame = cor.H_to_ref{matchedImgID} \ T;
                    global_frame = global_frame / global_frame(3,3);

                    world.frames_global(:, end+1) = ...
                        [
                        global_feature_id;
                        matchedImgID;
                        global_frame(1:2,3)
                        global_frame(1,1)
                        global_frame(2,1)
                        global_frame(1,2)
                        global_frame(2,2)                    
                        ];
                else
                    world.frames_global(:, end+1) = ...
                        [global_feature_id; matchedImgID; zeros(6,1)];
                end
                
            end
        end
        
        % Add unmatched features to map
        num_features = size(model.index.frames{matchedImgID}, 2);
        feature_vector = 1:num_features;
        unmatched_features = ~ismember(feature_vector, im2_feats);
        unmatched_features = feature_vector(unmatched_features);
        numUnmatchedFeats = length(unmatched_features);
        
        % Add new global feature for each unmatched local feature
        new_global_id = world.num_features + 1;
        world.num_features = new_global_id + numUnmatchedFeats - 1;
        new_global_ids = new_global_id : world.num_features;
        
        % Add link to local feature for each unmatched local feature
        for k = 1:numUnmatchedFeats
            world.feature_map(:, end+1) = ...
                [new_global_ids(k); matchedImgID; unmatched_features(k)];
            
            world.words_global(:, end+1) = ...
                [new_global_ids(k); matchedImgID; ...
                model.index.words{matchedImgID}(unmatched_features(k))];
            
            if ~isempty(cor.H_to_ref{matchedImgID})
                fndx = unmatched_features(k);

                % Transform local frames to global co-ordinate system
                % First obtain matrix transforming unit circle to
                % orientated elliptical frame in local co-ordinates
                T = [ ...
                    [model.index.frames{m}(3,fndx); model.index.frames{m}(4,fndx); 0] ...
                    [model.index.frames{m}(5,fndx); model.index.frames{m}(6,fndx); 0] ...
                    [model.index.frames{m}(1:2,fndx); 1] ...
                    ];
                % Then pre-multiply this by the transformation to global
                % co-ordinates to get the global feature ellipse
                global_frame = cor.H_to_ref{matchedImgID} \ T;
                global_frame = global_frame / global_frame(3,3);

                world.frames_global(:, end+1) = ...
                    [
                    new_global_ids(k);
                    matchedImgID;
                    global_frame(1:2,3)
                    global_frame(1,1)
                    global_frame(2,1)
                    global_frame(1,2)
                    global_frame(2,2)                    
                    ];
            else
                world.frames_global(:, end+1) = ...
                    [new_global_ids(k); matchedImgID; zeros(6,1)];
            end
            
        end
        
    end
    
end

% Need to work out the diff between frames and descriptors
% Might need to edit visualindex files to include exact descriptors in
% index
% At the moment using visual words instead of SIFT descriptors

end