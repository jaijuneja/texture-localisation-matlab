function [world, cor] = bundle_adjustment_local(world, cor, varargin)

opts.perspDistPenalty = 0;
opts.onlyOptimiseH = false;
opts.dH_thresh = 0.1;
opts.df_thresh = 1;
opts.imsToInclude = [];
opts = vl_argparse(opts, varargin);

% Note: we only optimise features that have been matched between multiple
% views, and then adjust unmatched features according to the new values of
% H_to_ref obtained

% Collect relevant features for optimisation
matched = (world.features_global(2,:) > 1) & world.features_mappable;
features_matched = world.features_global(:,matched);
indices_matched = world.feature_indices(:,matched);

% Collect relevant views for optimisation
ims_matched = find(cellfun(@(x)(~isempty(x)), cor.H_to_ref));
num_views = length(ims_matched);
num_feats_loc = nnz(indices_matched);

% NEW STUFF
if ~isempty(opts.imsToInclude)
    for i = 1:length(features_matched)
    matched_loc_features = indices_matched(:,i);
        for j = 1:nnz(matched_loc_features)
            imgid = world.frames_local(2, matched_loc_features(j));
            if ~ismember(imgid, opts.imsToInclude)
                indices_matched(j,i) = 0;
            end
        end
    end
    
    ims_matched = opts.imsToInclude;
    num_views = length(opts.imsToInclude);
    num_feats_loc = nnz(indices_matched);
end
%%%%

% Build a temporary cell array which contains the inverse values of
% cor.H_to_ref. We will perturb the inverted values during the optimisation
% and then re-invert them to obtain the new cor.H_to_ref
H_from_ref = cor.H_to_ref;
for i = 1:length(ims_matched)
    img = ims_matched(i);
    H_from_ref{img} = inv(cor.H_to_ref{img});
    H_from_ref{img} = H_from_ref{img} / H_from_ref{img}(3,3);
end

% Initialise the terms of the quadratic optimisation problem
A = zeros(num_feats_loc * 2, 1);
h = zeros(num_feats_loc * 2 + num_views * 8, 1);
G = zeros(num_feats_loc * 2 + num_views * 8);
if ~isequal(opts.perspDistPenalty, 0)
    A_penalty = zeros(num_feats_loc * 2, 1);
else
    A_penalty = [];
end

feat_counter = 0;
% Cycle through each global feature to build the optimisation parameters
for k = 1:length(features_matched)
    matched_loc_features = indices_matched(:,k);
    matched_loc_features = matched_loc_features(matched_loc_features ~= 0);
    for i = 1:length(matched_loc_features)
        feat_counter = feat_counter + 1;
        % Do a test here to check that the global features match
        frame_local = world.frames_local(:,matched_loc_features(i));
        
        % Determine position of feature and image in params A, h and G
        imgID = frame_local(2);
        imgLoc = find(ismember(ims_matched,imgID));

        img_ndx = (num_feats_loc * 2) + (imgLoc * 8) - 7;
        fLoc_ndx = feat_counter * 2 - 1;
        
        % Calculate parameters for quadprog
        % Recall the formula:
        % E(f, H) = ||?(f, H)||^2
        % Energy = a'*a + 2*a'*b'*delta + delta'*b*b'*delta
        % Quadprog takes form quadprog(G, h) for 1/2*x'*G*x + h'*x

        % Get the local optimisation parameters (a is the constant term, b
        % is the derivative term)
        [a, b] = get_optimisation_params(frame_local(3:4), ...
            features_matched(3:4,k), H_from_ref{imgID});
        
        % Plug the local values into the 'global' vector A
        A(fLoc_ndx:fLoc_ndx+1) = a;
        
        % Plug the local values into the 'global' vector h
        h_loc = 2 * a' * b'; % 10 x 1 vector
        h(fLoc_ndx:fLoc_ndx+1) = h(fLoc_ndx:fLoc_ndx+1) + h_loc(1:2)';
        h(img_ndx:img_ndx+7) = h(img_ndx:img_ndx+7) + h_loc(3:10)';

        % Plug the local values into the 'global' matrix G
        G_loc = (b * b'); % 10 x 10 matrix
        G(fLoc_ndx:fLoc_ndx+1,fLoc_ndx:fLoc_ndx+1) = ...
            G(fLoc_ndx:fLoc_ndx+1,fLoc_ndx:fLoc_ndx+1) + G_loc(1:2,1:2);
        G(fLoc_ndx:fLoc_ndx+1,img_ndx:img_ndx+7) = ...
            G(fLoc_ndx:fLoc_ndx+1,img_ndx:img_ndx+7) + G_loc(1:2,3:10);
        G(img_ndx:img_ndx+7,fLoc_ndx:fLoc_ndx+1) = ...
            G(img_ndx:img_ndx+7,fLoc_ndx:fLoc_ndx+1) + G_loc(3:10,1:2);
        G(img_ndx:img_ndx+7,img_ndx:img_ndx+7) = ...
        G(img_ndx:img_ndx+7,img_ndx:img_ndx+7) + G_loc(3:10,3:10);
    
        % We can introduce a penalty for perspective distortion, such that
        % Energy is increased if h31 and h32 deviate from zero
        if ~isequal(opts.perspDistPenalty, 0)
            [a31_penalty, h31_penalty, a32_penalty, h32_penalty] = ...
                constrain_perspective(H_from_ref{imgID}, opts.perspDistPenalty);
            h(img_ndx+2) = h(img_ndx+2) + h31_penalty;
            h(img_ndx+5) = h(img_ndx+5) + h32_penalty;
            G31_penalty = opts.perspDistPenalty;
            G32_penalty = opts.perspDistPenalty;
            G(img_ndx+2,img_ndx+2) = G(img_ndx+2,img_ndx+2) + G31_penalty;
            G(img_ndx+5,img_ndx+5) = G(img_ndx+5,img_ndx+5) + G32_penalty;
            % Test this in simple case first
            % Collect constant terms to check E before and after
            A_penalty(fLoc_ndx:fLoc_ndx+1) = [a31_penalty; a32_penalty];
        end
    end
end

A = [A; A_penalty];
E_before = A' * A;

G = 2 * G;

% Impose constraint that delta is small
H_ndxs = num_feats_loc*2+1:num_feats_loc*2 + num_views*8;
delta_upper = inf(size(h));
delta_upper(H_ndxs) = opts.dH_thresh;

F_ndxs = 1:num_feats_loc*2;
delta_upper(F_ndxs) = opts.df_thresh;

delta_lower = - delta_upper;

% If we're only optimising transformations (H), then remove all feature
% terms from optimisation parameters
if opts.onlyOptimiseH
    feat_ndxs = 1:num_feats_loc*2;
    G(feat_ndxs, :) = [];
    G(:, feat_ndxs) = [];
    h(feat_ndxs) = [];
    delta_upper = [];
    delta_lower = [];
    num_feats_loc = 0;
end

options = optimoptions('quadprog');
options.Algorithm = 'trust-region-reflective';
options.Display = 'off';
delta = quadprog(G, h, [], [], [], [], delta_lower, delta_upper);
% delta = fmincon(FUN(delta),X0,A,B,Aeq,Beq,LB,UB,NONLCON)

% delta(num_feats_loc * 2 + 1:end) = 0; % Ignore transformations (for testing purposes)
% delta(1:num_feats_loc * 2) = 0; % Ignore features (for testing purposes)

E_after = A'*A + h'*delta + 0.5*delta'*G*delta;

%%%%%%%%%%%%%%%%%%%%% START ENERGY CALCULATION TEST %%%%%%%%%%%%%%%%%%%%%%%
% Calculate energy using naive approach and compare with result above
doEnergyCalcTest = 0;
if doEnergyCalcTest
    E_current = 0;
    h31_ndx = 5;
    h32_ndx = 8;
    feat_counter = 0;
    for k = 1:length(features_matched)
        matched_loc_features = indices_matched(:,k);
        matched_loc_features = matched_loc_features(matched_loc_features ~= 0);
        for i = 1:length(matched_loc_features)
            feat_counter = feat_counter + 1;
            fLoc_ndx = feat_counter * 2 - 1;
            
            frame_local = world.frames_local(:,matched_loc_features(i));
            
            % Determine index of transformation in global optimisation params
            imgID = frame_local(2);
            imgLoc = find(ismember(ims_matched,imgID));
            img_ndx = (num_feats_loc * 2) + (imgLoc * 8) - 7;
            
            % Get the local optimisation parameters
            [a, b] = get_optimisation_params(frame_local(3:4), ...
                features_matched(3:4,k), H_from_ref{imgID});
            
            if opts.onlyOptimiseH
                b(1:2,:) = [];
                fLoc_ndx = [];
                h31_ndx = 3;
                h32_ndx = 6;
            end
            
            % Collect the local delta values into a new delta vector
            delta_new = [delta(fLoc_ndx:fLoc_ndx+1); delta(img_ndx:img_ndx+7)];
            
            % Incrementally calculate energy at each iteration
            E_current = E_current + double(norm(a + b' * delta_new))^2;
            
            % Add energy due to perspective distortion penalty
            E_current = E_current + opts.perspDistPenalty * ...
                ( (H_from_ref{imgID}(3,1) + delta_new(h31_ndx))^2 + ...
                (H_from_ref{imgID}(3,2) + delta_new(h32_ndx))^2 );
        end
    end
    % Check that naive calculation yields same result
    assert(abs((E_current - E_after)/E_current) < 1e-3) ;
end
%%%%%%%%%%%%%%%%%%%%%%% END ENERGY CALCULATION TEST %%%%%%%%%%%%%%%%%%%%%%%

% Update features
feature_ndxs = nonzeros(indices_matched);
if ~opts.onlyOptimiseH
    f_delta = reshape(delta(1:num_feats_loc*2), 2, num_feats_loc);
    world.frames_local(3:4,feature_ndxs) = world.frames_local(3:4,feature_ndxs) + f_delta;
    % Do we sort out the elliptical bit? If we know the new H, we can
    % calculate a new ellipse from the change in H, perhaps?
end

% Update global transformations
for i = 1:num_views
    imgID = ims_matched(i);
    img_ndx = num_feats_loc*2+i*8-7;
    for j = 1:8
        H_from_ref{imgID}(j) = H_from_ref{imgID}(j) + delta(img_ndx+j-1);
    end
    cor.H_to_ref{imgID} = inv(H_from_ref{imgID});
    cor.H_to_ref{imgID} = cor.H_to_ref{imgID}/cor.H_to_ref{imgID}(3,3);
end

% We need a version of build_world that only updates single-view features
% Maybe a function update_world that asks what information is used for
% update - e.g. use cor.H_to_ref, ignore matched feats. Use cor.H, ignore
% no feats.

fprintf(['Iteration of bundle adjustment completed: \n ' ...
    'Energy before: %d \n Energy after: %d \n'], E_before, E_after);

%%%%%%%%%%%%%%%%% START LINEARIZATION ERROR CALCULATION %%%%%%%%%%%%%%%%%%%
% Calculate energy using naive approach and compare with result above
doLinearErrorTest = 0;
if doLinearErrorTest
    H_from_ref = cor.H_to_ref;
    for i = 1:length(ims_matched)
        img = ims_matched(i);
        H_from_ref{img} = inv(cor.H_to_ref{img});
        H_from_ref{img} = H_from_ref{img} / H_from_ref{img}(3,3);
    end
    E_after_true = 0;
    for k = 1:length(features_matched)
        matched_loc_features = indices_matched(:,k);
        matched_loc_features = matched_loc_features(matched_loc_features ~= 0);
        for i = 1:length(matched_loc_features)
            
            frame_local = world.frames_local(:,matched_loc_features(i));
            
            % Determine index of transformation in global optimisation params
            imgID = frame_local(2);
            
            % Get the local optimisation parameters
            [a, ~] = get_optimisation_params(frame_local(3:4), ...
                features_matched(3:4,k), H_from_ref{imgID});
            
            % Incrementally calculate energy at each iteration
            E_after_true = E_after_true + double(norm(a))^2;
            
            E_after_true = E_after_true + opts.perspDistPenalty * ...
                ( H_from_ref{imgID}(3,1)^2 + H_from_ref{imgID}(3,2)^2 );
        end
    end
    lin_error_pc = round(abs((E_after_true - E_after)/E_after_true)*100);
    % Check that naive calculation yields same result
    fprintf(['Actual energy after optimisation is %d \n (%d%% error ' ...
        'from linearised energy prediction) \n\n'], E_after_true, lin_error_pc);
end
%%%%%%%%%%%%%%%%%% END LINEARIZATION ERROR CALCULATION %%%%%%%%%%%%%%%%%%%%

end