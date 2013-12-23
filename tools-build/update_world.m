function world = update_world(world, cor)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 20/12/2013
% -------------------------------------------------------------------------
%
% UPDATE_WORLD
% world = build_world(model, cor)
% 
% Given new transformations cor.H_to_world, update the positions of all
% unmatched features.
%
% Inputs:
%   - world:    World structure containing global features. Type 'help 
%               build_world' for more info
%   - cor:      Correspondence structure containing links between different
%               images (graph representation). Type
%               'help build_correspondence' for more info
%
% Outputs:
%   - world

% Get the ids of all unmatched global features
unmatched_global = find(world.features_global(2,:) == 1);
numUnmatchedFeats = length(unmatched_global);

% Then get their corresponding local feature ids
unmatched_local = world.feature_indices(1, unmatched_global);

% Given local and global feature ids, update the position of unmatched
% global feats
for i = 1:numUnmatchedFeats
    global_id = unmatched_global(i);
    local_id = unmatched_local(i);
    img_id = world.frames_local(2, local_id);
    
    local_frame = world.frames_local(3:8, local_id);
    
    % Update global feature frame
    global_frame = transform_frames(local_frame, cor.H_to_ref{img_id});
    world.features_global(3:4, global_id) = global_frame(1:2);
end

end