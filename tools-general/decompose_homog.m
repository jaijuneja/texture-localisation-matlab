function [R, t] = decompose_homog(H, K)

r1r2t = K \ H;

r1 = r1r2t(:,1);
scale1 = norm(r1);
% r1 = r1 / scale1;
r2 = r1r2t(:,2);
scale2 = norm(r2);
% r2 = r2/scale2;
% r2 = r2/norm(r2);
r3 = cross(r1, r2);
% r3 = r3/norm(r3);
% r1 = cross(r2, r3);

R = [r1 r2 r3];
% Improvement to commented out method above - see planar-struct3.pdf
[U, ~, V] = svd(R);
R = U * V';

% Want to scale t so that it represents r1r2t as closely as possible
scalet = mean([scale1 scale2]);
% Alternative scaling
% scalet2 = norm(r1r2t(:,1:2))/norm(R(:,1:2));
t = r1r2t(:,3);
t = t/scalet;

end