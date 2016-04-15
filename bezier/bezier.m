function out = bezier(coords, t)
%BEZIER Calculates cubic bezier curve based on control points
%
% out = bezier(coords)
% out = bezier(coords, t)
%
% This function calculates a cubic bezier curve based on the input control
% points.  If more than 4 points are given, the points are divided into
% sets of 4 (overlapping at endpoints), with the resulting curves patched
% together for output.
%
% Input variables:
%
%   coords: n x m array of control point coordinates, where n is the number
%           of points and must be from the set 4+3x, where x is an integer,
%           and m is the number of dimensions needed to describe the point
%           locations.
%
%   t:      vector of values between 0 and 1 used to calculate each bezier
%           curve.  If not included, default is t = linspace(0,1,101)
%
% Output variables:
%
%   out:    length(t) x m array of output coordinates defining the bezier
%           curve 

% Copyright 2007 Kelly Kearney


%-------------------------
% Check input
%-------------------------

narginchk(1,2);

npoints = size(coords,1);

if mod(npoints - 4, 3) ~= 0
    error('Number of points must 4 + 3n');
end

if nargin == 1
    t = linspace(0,1,101);
end

%-------------------------
% Calculate beziers
%-------------------------

ncurve = (npoints - 4)/3 + 1;

out = [];


for icurve = 1:ncurve

    istart = (icurve-1)*3 + 1;
    
    % Equation of Bezier Curve, utilizes Horner's rule for efficient 
    % computation.
    % Q(t)=(-P0 + 3*(P1-P2) + P3)*t^3 + 3*(P0-2*P1+P2)*t^2 + ...
    %      3*(P1-P0)*t + Px0
    
    p0 = coords(istart,:);
    p1 = coords(istart+1,:);
    p2 = coords(istart+2,:);
    p3 = coords(istart+3,:);
    
    c3 = -p0 + 3*(p1-p2) + p3;
    c2 = 3*(p0 - (2*p1)+p2); 
    c1 = 3*(p1 - p0);
    c0 = p0;
    for k=1:length(t)
        q(k,:)=((c3*t(k)+c2)*t(k)+c1)*t(k) + c0;    
    end
    
    % Remove duplicate points and add to output
    
    out = [out; q];

end

% Remove duplicate points

% isdup = ~any(diff(out),2);
% out(isdup,:) = [];


