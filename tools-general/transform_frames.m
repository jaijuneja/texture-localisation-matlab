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
    new_frame = H * T;
    new_frame = new_frame / new_frame(3,3);
    
    if ~(isequal(new_frame(3,1), 0) && isequal(new_frame(3,2), 0))
        % Frame is not affine - need to force affinity
        t_loc = frames(1:2, i); % new_frame(1:2, 3);
        t_glob = new_frame(1:2, 3);
        A = get_jacobian(new_frame, t_loc);
        new_frame = [A t_glob; 0 0 1];
    end
        
    new_frames(:,i) = [ new_frame(1:2,3)
                        new_frame(1,1)
                        new_frame(2,1)
                        new_frame(1,2)
                        new_frame(2,2) ];
end

function jac = get_jacobian(H, t)
x = t(1);
y = t(2);

jac = [ ...
    [ H(1,1)/(H(3,3) + H(3,1)*x + H(3,2)*y) - ...
    (H(3,1)*(H(1,3) + H(1,1)*x + H(1,2)*y))/(H(3,3) + H(3,1)*x + H(3,2)*y)^2 ; ...
    H(2,1)/(H(3,3) + H(3,1)*x + H(3,2)*y) - ...
    (H(3,1)*(H(2,3) + H(2,1)*x + H(2,2)*y))/(H(3,3) + H(3,1)*x + H(3,2)*y)^2 ] ...
    [ H(1,2)/(H(3,3) + H(3,1)*x + H(3,2)*y) - ...
    (H(3,2)*(H(1,3) + H(1,1)*x + H(1,2)*y))/(H(3,3) + H(3,1)*x + H(3,2)*y)^2 ; ...
    H(2,2)/(H(3,3) + H(3,1)*x + H(3,2)*y) - ...
    (H(3,2)*(H(2,3) + H(2,1)*x + H(2,2)*y))/(H(3,3) + H(3,1)*x + H(3,2)*y)^2 ] ];
 
% function jac = get_jacobian_sym
% syms x y xnew ynew h11 h12 h13 h21 h22 h23 h31 h32 h33;
% 
% F(x, y) = [(h11 * x + h12 * y + h13) / (h31 * x + h32 * y + h33); ...
%     (h21 * x + h22 * y + h23) / (h31 * x + h32 * y + h33)];
% 
% dF_dx = diff(F, x);
% dF_dy = diff(F, y);