function new_frames = transform_frames(frames, H)
% Jai Juneja, www.jaijuneja.com
% University of Oxford
% 08/12/2013
% -------------------------------------------------------------------------
%
% TRANSFORM_FRAMES
% new_frames = transform_frames(frames, H)
%
% Transforms oriented elliptical SIFT frames using the transformation
% matrix H.
%
% Inputs:
%   - frames:   6xn matrix for n frames, where frames(1:2,:) describes the
%               (x,y) co-ordinates of the features; frames(3:6) describes
%               elements A(1,1), A(2,1), A(1,2) and A(2,2) respectively of
%               the matrix that transforms the unit circle to the SIFT
%               ellipse in the local frame
%   - H:        3x3 transformation matrix
%
% Outputs:
%   - new_frames:   Matrix of transformed frames; has the same dimensions
%                   as input argument 'frames'

new_frames = zeros(size(frames));
for i = 1:size(frames,2)
    % First obtain matrix transforming unit circle to orientated elliptical
    % frame in local co-ordinates
    T = [ ...
        [frames(3,i); frames(4,i); 0] ...
        [frames(5,i); frames(6,i); 0] ...
        [frames(1:2,i); 1] ];

    % Then pre-multiply this by the transformation H to get the transformed
    % feature ellipse
    global_frame = H * T;
    global_frame = global_frame / global_frame(3,3);
    
    new_frames(:,i) = [ global_frame(1:2,3)
                        global_frame(1,1)
                        global_frame(2,1)
                        global_frame(1,2)
                        global_frame(2,2) ];
end

end