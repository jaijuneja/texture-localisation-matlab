function world = update_features(world, cor, model, ...
    fGlobalID, imgID, fLocalID)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/11/2013
% -------------------------------------------------------------------------
% 
% UPDATE_FEATURES
% world = update_features(world, cor, model, fGlobalID, imgID, fLocalID)
%
% Given a global feature fGlobalID and a local feature fLocalID in an image
% imgID to be added to the world, update the feature map, words and frames
% in the the world structure.
%
% Inputs:
%   - world:        World structure containing global features. Type 'help 
%                   build_world' for more info
%   - cor:          Correspondence structure containing links between
%                   different images. Type 'help build_correspondence' for
%                   more info
%   - model:        Index of images from visualindex. Type 'help
%                   visualindex_build' for more info
%   - fGlobalID:    Global ID of feature to be added (could be a new ID)
%   - imgID:        Image ID from which the feature is being added to the
%                   global map
%   - fLocalID:     Local ID of feature to be added from image imgID
%
% Outputs:
%   - world

% Insert the feature into the global map
world.feature_map(:, end+1) = ...
    [fGlobalID; imgID; fLocalID];

% Increment number of features in global map
world.num_features = world.num_features + 1;

% Also update words in same way
world.words_global(:, end+1) = ...
    [fGlobalID; imgID; ...
    model.index.words{imgID}(fLocalID)];

if ~isempty(cor.H_to_ref{imgID})
    % Transform local frames to global co-ordinate system
    % First obtain matrix transforming unit circle to
    % orientated elliptical frame in local co-ordinates
    T = [ ...
        [model.index.frames{imgID}(3,fLocalID); ...
        model.index.frames{imgID}(4,fLocalID); 0] ...
        [model.index.frames{imgID}(5,fLocalID); ...
        model.index.frames{imgID}(6,fLocalID); 0] ...
        [model.index.frames{imgID}(1:2,fLocalID); 1] ...
        ];
    % Then pre-multiply this by the transformation to global
    % co-ordinates to get the global feature ellipse
    global_frame = cor.H_to_ref{imgID} \ T;
    global_frame = global_frame / global_frame(3,3);

    world.frames_global(:, end+1) = ...
        [   fGlobalID
            imgID
            global_frame(1:2,3)
            global_frame(1,1)
            global_frame(2,1)
            global_frame(1,2)
            global_frame(2,2)   ];
    
    % The feature's global co-ordinates are known - it can be mapped
    world.features_mappable(end+1) = true;
    
else
    world.frames_global(:, end+1) = ...
        [fGlobalID; imgID; nan(6,1)];
    
    % The feature's global co-ordinates are unknown - it cannot be mapped
    world.features_mappable(end+1) = false;
end

end