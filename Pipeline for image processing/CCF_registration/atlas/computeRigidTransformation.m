function [ transformation ] = computeRigidTransformation( points1, points2 )
%COMPUTERIGIDTRANSFORMATION Computers a rigid transformation from points1
%to points2
%   This functions assumes that all points are inliers and uses SVD to
%   compute the best transformation. 
%  modified by YZ, for matlab coordinate
    nPoints = size(points1, 1);
    dimension = size(points1, 2);
    
    % previous: center
%     centroid1 = sum(points1, 1)./nPoints;
%     centroid2 = sum(points2, 1)./nPoints;

    % current: 0, 0
    centroid1 = [1, 1]; % centroids for x
    centroid2 = [1, 1]; % centroids for y

    
    
    centered1 = points1 - repmat(centroid1, nPoints, 1); % decentered?
    centered2 = points2 - repmat(centroid2, nPoints, 1); % decentered?
    
    W = eye(nPoints, nPoints);
    
    S = centered1' * W * centered2;
    
    [U, Sigma, V] = svd(S);
    
    M = eye(dimension, dimension);
    M(dimension, dimension) = det(V*U');
    
    R = V * M * U';
    t = centroid2' - R*centroid1';
    
    % build the transformation matrix
    transformation = eye(3,3);
    transformation(1:2, 1:2) = R; % rotation
    transformation(1:2, 3) = t;
%     transformation(3, 2 : -1 : 1) = t; % translation
end
