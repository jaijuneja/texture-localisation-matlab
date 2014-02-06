function x = transform_world(x, T)

if isequal(size(T), [1 1])
    % T is a scalar that represents a scale factor
    scaling_mat = eye(3);
    scaling_mat(3,3) = 1/T;
    T = sparse(scaling_mat);
end

if isstruct(x)
    if isfield(x, 'H_to_world')
        % x is the correspondence structure
        mappable = find(cellfun(@(x)(~isempty(x)), x.H_to_world));
        for i = mappable
            x.H_to_world{i} = T * x.H_to_world{i};
            x.H_to_world{i} = x.H_to_world{i} * x.H_to_world{i}(3,3);
        end
        x.H_world_toref = x.H_world_toref / T;
        x.H_world_toref = x.H_world_toref / x.H_world_toref(3,3);
        
    elseif isfield(x, 'features_global')
        % x is the world structure
        x.features_global(3:4,:) = transform_points( ...
            x.features_global(3:4,:), T);
    end
else
    % x is just a transformation matrix
    x = T * x;
    x = x / x(3,3);
end

end