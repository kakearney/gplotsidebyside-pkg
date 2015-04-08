function [h, A] = gplotsidebyside(adj, varargin)
%GPLOTSIDEBYSIDE Graph with edges lined up side-by-side around a circle
%
% [h, A] = gplotsidebyside(adj, p1, v1, ...)
%
% Input variables:
%
%   adj:        n x n x m adjacency matrix; 3rd dimension allows for
%               multiple edges between the same set of nodes 
%
% Optional input variables:
%
%   nodelabel:  n x 1 cell array of strings {'1';'2';'3';...}
%
%   rout:       radius of outer circle [1.0]
%
%   rin:        radius of inner circle (used as control points for bezier
%               curves) [0.25]
%
%   rlabel:     1 x 2 array, inner and outer radius of node separation 
%               ticks [1.0 1.1]
%
%   rnode:      1 x 2 array, inner and outer radius of nodes
%
%   nb:         number of points along each edge curve [20]
%
%   nval:       nnode x 1 array, values corresponding to nodes, to be
%               scaled in same way as edges
%
% Output variables:
%
%   h:          structure of graphics handles
%
%   A:          structure of data


% Copyright 2014 Kelly Kearney


%-----------------
% Options
%-----------------

nnode = size(adj,1);

Opt.nodelabel = strtrim(cellstr(num2str((1:nnode)')));
Opt.rout = 1;
Opt.rin = 0.25;
Opt.rlabel = [1.0 1.1];
Opt.axis = gca;
Opt.nval = [];
Opt.rnode = [1.025 1.075];
Opt.nb = 20;

Opt = parsepv(Opt, varargin);

%-----------------
% Calculations
%-----------------

% Links

adj = full(adj);

idx = find(adj);
[ii,jj,kk] = ind2sub(size(adj), idx);

links = sortrows([ii jj kk adj(adj~=0)]);

% For each node, links are arranged with all inputs, then all outputs

nlink = nnz(adj);

tmp = cat(2, [ii jj kk; jj ii kk], [zeros(nlink,1); ones(nlink,1)], repmat(adj(adj~=0),2,1));
[srt, isrt] = sortrows(tmp, [1 4 2]);

lims = minmax(srt(:,end), 'expand', 0.01);
if lims(1) == lims(2)
    relwidth = ones(size(srt(:,end)));
else
    relwidth = interp1(lims, [0 1], srt(:,end));
end

relsum = sum(relwidth);
theta = ([0; cumsum(relwidth)]./relsum) * 2 * pi;
% theta = [0; cumsum(relwidth)];
% theta = theta./(max(theta)) * 2*pi;
dtheta = diff(theta);
thetamid = (theta(1:end-1)+theta(2:end))./2;

if ~isempty(Opt.nval)
    relwidthnode = interp1(lims, [0 1], Opt.nval, 'linear', 'extrap');
    dthetanode = (relwidthnode./relsum) * 2 * pi;
    dthetanode = reshape(dthetanode,1,[]);
end

% Disentangle which is which new subnode corresponds to the source and sink
% of each link

srt2 = srt;
srt3 = srt;
isr = srt(:,4) == 1;
srt2(isr,:) = NaN;
srt3(~isr,:) = NaN;

[tf, src] = ismember(links, srt2(:,[1 2 3 5]), 'rows'); 
[tf, snk] = ismember(links, srt3(:,[2 1 3 5]), 'rows'); 

% Calculate edge curve coordinates

xout = cos(theta) * Opt.rout;
yout = sin(theta) * Opt.rout;

xin = cos(theta) * Opt.rin;
yin = sin(theta) * Opt.rin;

adjtmp = sparse(src, snk, links(:,4), nlink*2, nlink*2);
xmid = cos(thetamid) * Opt.rout;
ymid = sin(thetamid) * Opt.rout;

tb = linspace(0,1,Opt.nb);

[xedge,yedge] = deal(zeros(Opt.nb*2, nlink));

for il = 1:nlink
    
    x1 = [xout(src(il)) xin(src(il)) xin(snk(il)+1) xout(snk(il)+1)];
    x2 = [xout(src(il)+1) xin(src(il)+1) xin(snk(il)) xout(snk(il))];
    
    y1 = [yout(src(il)) yin(src(il)) yin(snk(il)+1) yout(snk(il)+1)];
    y2 = [yout(src(il)+1) yin(src(il)+1) yin(snk(il)) yout(snk(il))];
    
    xy1 = bezier([x1' y1'], tb);
    xy2 = bezier([x2' y2'], tb);
    
    xypoint = (xy1(end,:) + xy2(end,:))./2;
    
    xy = [xy1(1:end-1,:); xypoint; xy2(end-1:-1:1,:); xy1(1,:)];
%     xy = [xy1; flipud(xy2); xy1(1,:)];
    
    xedge(:,il) = xy(:,1);
    yedge(:,il) = xy(:,2);
    
end

% Coordinates of lines separating each group of fluxes

[node,tmp] = aggregate(srt(:,1), theta(2:end), @max);
thsep = nan(nnode,1);
thsep(node) = cell2mat(tmp);

thmid = ([0; thsep(1:end-1)]+thsep)./2;

xsep = cos(thsep) * Opt.rlabel;
ysep = sin(thsep) * Opt.rlabel;

xtext = cos(thmid) * mean(Opt.rlabel);
ytext = sin(thmid) * mean(Opt.rlabel);

if ~isempty(Opt.nval)
    thnode = [thmid thmid] + dthetanode' * [-0.5 0.5];
    
    x1 = cos(thnode) * Opt.rnode(1);
    x2 = cos(thnode) * Opt.rnode(2);
    y1 = sin(thnode) * Opt.rnode(1);
    y2 = sin(thnode) * Opt.rnode(2);

    xnode = [x1(:,1:2) x2(:,[2 1]) x1(:,1)];
    ynode = [y1(:,1:2) y2(:,[2 1]) y1(:,1)];
    
end

%-----------------
% Plot
%-----------------

% Plot edges, shaded by relative value

axes(Opt.axis);
hold(Opt.axis, 'on');

h.edge = patch(xedge, yedge, srt(src,end)');
axis equal;
set(Opt.axis, 'visible', 'off');

if ~isempty(Opt.nval)
    h.node = patch(xnode', ynode', Opt.nval);
end

h.sep = plot(xsep', ysep', 'k');
h.txt = text(xtext, ytext, Opt.nodelabel, 'horiz', 'center');

setappdata(h.edge, 'index', srt(src,1:3));

A.lims = lims;
A.thmid = thmid;
A.edgeindex = srt(src,1:3);








