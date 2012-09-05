function handles = plotmod(tr,varargin)
%PLOT draws a phylogenetic tree.
%
%   PLOT(TREE) draws a static plot of a phylogenetic tree object. PLOT
%   allows further customization of a rendered tree either by using
%   parameter-value pairs or by manually modifying the HG properties of the 
%   graph elements. For an interactive exploration tool use instead the
%   VIEW method.  
%
%   PLOT(TREE,ACTIVEBRANCHES) selects a subset of the branches in the
%   tree for the plot, all other branches are not drawn. ACTIVEBRANCHES is 
%   a logical array of size [numBranches x 1] indicating the active
%   branches. 
%
%   H = PLOT(...) returns a structure with handles to the graph elements.
%   Handles to graph elements are also stored in the 'UserData' figure
%   field. 
%
%   PLOT(...,'TYPE',TYPE) selects the method to render the phylogenetic
%   tree. Options are: 'square' (default), 'angular', 'radial',
%   'equalangle', and 'equaldaylight'.
%
%   PLOT(...,'ORIENTATION',POS) sets the position of the root node of the
%   phylogenetic tree relative to the axes. Options are: 'top', 'bottom',
%   'left' (default), and, 'right'. Orientation parameter is  valid only
%   when TYPE is either 'square' or 'angular'.
%
%   PLOT(...,'ROTATION',ANGLE) sets the angle (in degrees) that rotates the
%   phylogenetic tree. ANGLE is a scalar between 0 and 360 and defaults to
%   0. Rotate parameter is only valid when TYPE is either 'radial',
%   'equalangle', or 'equaldaylight'.
%
%   PLOT(...,'BRANCHLABELS',TF) displays the branch labels when TF is true.
%   Default is false. Branch labels are placed next to the branch node.
%
%   PLOT(...,'LEAFLABELS',TF) displays the leaf labels when TF is true. 
%   Default is false when TYPE is either 'square' or 'angular', otherwise
%   default is true. Leaf labels are placed next to the leaf nodes.
%
%   PLOT(...,'TERMINALLABELS',TF) hides the terminal labels when TF is
%   false. Default is true. Terminal labels are placed at the axis as
%   tick labels. Terminal labels can only be drawn when TYPE is either
%   'square' or 'angular'.
%
%   PLOT(...,'LLROTATION',TF) rotates leaf labels so that the text is
%   aligned to the root when TF is true. Default is false. Leaf labels can
%   only be rotated when TYPE is either 'radial', 'equalangle', or
%   'equaldaylight'. 
%
%   Example:
%
%       tr = phytreeread('pf00002.tree')
%       h = plot(tr,'type','radial')
%
%       % Graph element properties can be modified as follows:
%       set(h.branchNodeLabels,'FontSize',7,'Color',[0 0 .5])
%
%   See also PHYTREE, PHYTREE/VIEW, PHYTREEREAD, PHYTREETOOL, SEQNEIGHJOIN,
%   SEQLINKAGE. 

% Copyright 2003-2009 The MathWorks, Inc.
% $Revision: 1.1.6.13 $ $Author: batserve $ $Date: 2009/05/07 18:15:51 $

if numel(tr)~=1
    error('Bioinfo:phytree:plot:NoMultielementArrays',...
        'Phylogenetic tree must be an 1-by-1 object.');
end

% set defaults
dispBranchLabels = false;
dispLeafLabels = NaN;
dispTerminalLabels = true;
renderType = 'square';
circularType = false;
orientation = 'left';
rotation = 0;
rotateLeafLabels = false;
% Default fontsize for labels 
fontSizeLeaf = 8;
fontSizeBranch = 8;
fontSizeTerminal = 9;

leafDisplayedAsBranch = [];

if nargin>1 && islogical(varargin{1})
    oldLeafNames = tr.names(1:get(tr,'numLeaves'));
    tr = prune(tr,~varargin{1},'mode','exclusive');
    prunedNames = get(tr,'LeafNames');
    ws = warning('off','Bioinfo:seqmatch:StringNotFound');
    order = seqmatch(prunedNames,oldLeafNames);
    warning(ws)
    leafDisplayedAsBranch = find(order==0)';
    argStart = 2;
else
    argStart = 1;
end

tr = struct(tr);
tr.numBranches = size(tr.tree,1);
tr.numLeaves = tr.numBranches + 1;
tr.numLabels = tr.numBranches + tr.numLeaves;

if nargin - argStart > 0
    if rem(nargin - argStart,2) == 1
        error('Bioinfo:phytree:plot:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'type','orientation','rotation',...
        'branchlabels','leaflabels','terminallabels',...
        'llrotation'};
    for j = argStart:2:nargin-argStart
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:phytree:plot:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:phytree:plot:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % type
                    oktypes={'square','angular','radial',...
                             'equalangle','equaldaylight'};
                    l = strmatch(lower(pval),oktypes); 
                    if isempty(l)
                        error('Bioinfo:phytree:plot:UnknownTypeName',...
                            'Unknown option for %s.',upper(okargs{k}));
                    elseif length(l)>1
                        error('Bioinfo:phytree:plot:AmbiguousTypeName',...
                            'Ambiguous option for %s.',upper(oktypes{l}));
                    else
                        renderType = oktypes{l};
                        if l>2
                            circularType = true;
                        end
                    end
                case 2 % orientation
                    oktypes={'left','right','top','bottom'};
                    l = strmatch(lower(pval),oktypes); 
                    if isempty(l)
                        error('Bioinfo:phytree:plot:UnknownOrientation',...
                            'Unknown option for %s.',upper(okargs{k}));
                    else
                        orientation = oktypes{l};
                    end
                case 3 % rotation
                    if isreal(pval(1))
                        rotation = double(pval(1))*pi/180;
                    else
                        error('Bioinfo:phytree:plot:NotValidType',...
                            'ROTATION must be numeric and real');
                    end
                case 4 % branch labels
                    dispBranchLabels = opttf(pval);
                case 5 % leaf labels
                    dispLeafLabels = opttf(pval);
                case 6 % terminal labels
                    dispTerminalLabels = opttf(pval);
                case 7 % rotate leaf labels
                    rotateLeafLabels = opttf(pval);
            end
        end
    end
end

% set dependent default for dispLeafLabels
if isnan(dispLeafLabels)
    if circularType
        dispLeafLabels = true;
    else
        dispLeafLabels = false;
    end
end
% terminal labels are not allowed on circular layouts
if circularType
    dispTerminalLabels = false;
end

% empty names default to generic names
for ind = 1:tr.numLabels;
    if isempty(tr.names{ind})
        if ind > tr.numLeaves
            tr.names{ind} = ['Branch ' num2str(ind-tr.numLeaves)];
        else
            tr.names{ind} = ['Leaf ' num2str(ind)];
        end
    end
end

% some helper variables
nodeIndex   = 1:tr.numLabels;
leafIndex   = 1:tr.numLeaves;
branchIndex = tr.numLeaves+1:tr.numLabels;
prunedLeafIndex = setdiff(leafIndex,leafDisplayedAsBranch);
prunedBranchIndex = [branchIndex leafDisplayedAsBranch];

% obtain parents for every node
tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

% propagate last leaf
tr.lastleaf = (1:tr.numLabels)';

% find x-y coordinates of nodes in tree according with the layout type
tr.terminalNodes = tr.lastleaf([true;diff(tr.lastleaf(1:tr.numLeaves))~=0]);
tr.y=zeros(tr.numLabels,1);
tr.y(tr.terminalNodes)=1:length(tr.terminalNodes);
switch renderType
    case {'square','angular','radial'}
        % find x coordinates of branches:
        tr.x = tr.dist;
        for ind = tr.numBranches:-1:1
            tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
            
        end
        tr.x=tr.x+(1-max(tr.x));
        % find y coordinates of terminal nodes:
        tr.terminalNodes = tr.lastleaf([true;diff(tr.lastleaf(1:tr.numLeaves))~=0]);
        tr.y=zeros(tr.numLabels,1);
        tr.y(tr.terminalNodes)=1:length(tr.terminalNodes);
        % propagate proper y coordinates to the branch nodes:
        if strcmp(renderType,'angular')
            for ind = 1:tr.numBranches
                if tr.x(tr.tree(ind,1))/tr.x(tr.tree(ind,2))>3
                    tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,1));
                elseif tr.x(tr.tree(ind,2))/tr.x(tr.tree(ind,1))>3
                    tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,2));
                else
                    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
                end
            end
        else % 'square','radial'
            for ind = 1:tr.numBranches
                tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
            end
        end
    case{'equalangle', 'equaldaylight'}
        tr.nchild = sum(tr.tree <= tr.numLeaves, 2);
        for ind = 2:tr.numBranches
            treerow = tr.tree(ind,:);
            tr.nchild(ind) = tr.nchild(ind) + sum(tr.nchild(treerow(treerow > tr.numLeaves) - tr.numLeaves));
        end
        tr.nchild = [ones(tr.numLeaves, 1); tr.nchild];
        tr.wedge = tr.nchild*2*pi/tr.numLeaves;
        
        tr.ang = zeros(tr.numLabels, 1);
        tr.x = zeros(tr.numLabels,1);
        tr.y = zeros(tr.numLabels,1);
        for i = tr.numBranches:-1:1
            parNode = i + tr.numLeaves;
            rChild = tr.tree(i, 1);
            lChild = tr.tree(i, 2);
            
            tr.ang(rChild) = tr.ang(parNode) + tr.wedge(lChild)/2;
            tr.ang(lChild) = tr.ang(parNode) - tr.wedge(rChild)/2;
            
            tr.x(rChild) = tr.x(parNode) + tr.dist(rChild)*cos(tr.ang(rChild));
            tr.y(rChild) = tr.y(parNode) + tr.dist(rChild)*sin(tr.ang(rChild));
            
            tr.x(lChild) = tr.x(parNode) + tr.dist(lChild)*cos(tr.ang(lChild));
            tr.y(lChild) = tr.y(parNode) + tr.dist(lChild)*sin(tr.ang(lChild));
        end
        
        tr.tang = tr.ang;
        
        if strcmp(renderType, 'equaldaylight')
            for selectedNode = branchIndex(end-1:-1:1)
                
                subtree1 = getChildNodes(tr.tree(selectedNode - tr.numLeaves, 1), tr);
                subtree2 = getChildNodes(tr.tree(selectedNode - tr.numLeaves, 2), tr);
                subtree3 = nodeIndex(~ismember(nodeIndex, [subtree1 subtree2 selectedNode]));
                
                subtree1Leaf = subtree1(subtree1 <= tr.numLeaves);
                subtree2Leaf = subtree2(subtree2 <= tr.numLeaves);
                subtree3Leaf = subtree3(subtree3 <= tr.numLeaves);
                
                xdiff = tr.x(tr.tree(selectedNode - tr.numLeaves, :)) - tr.x(selectedNode);
                ydiff = tr.y(tr.tree(selectedNode - tr.numLeaves, :)) - tr.y(selectedNode);
                
                refangles(1:2) = atan2(ydiff, xdiff);
                refangles(3) = atan2(tr.y(tr.par(selectedNode)) - tr.y(selectedNode), ...
                    tr.x(tr.par(selectedNode)) - tr.x(selectedNode));
                
                st1angles = atan2((tr.y(subtree1Leaf) - tr.y(selectedNode)), (tr.x(subtree1Leaf) - tr.x(selectedNode))) - refangles(1);
                st2angles = atan2((tr.y(subtree2Leaf) - tr.y(selectedNode)), (tr.x(subtree2Leaf) - tr.x(selectedNode))) - refangles(2);
                st3angles = atan2((tr.y(subtree3Leaf) - tr.y(selectedNode)), (tr.x(subtree3Leaf) - tr.x(selectedNode))) - refangles(3);
                
                b = (st1angles < -pi) | (st1angles > pi);
                st1angles(b) = st1angles(b) - floor(st1angles(b)./pi)*2*pi;
                b = (st2angles < -pi) | (st2angles > pi);
                st2angles(b) = st2angles(b) - floor(st2angles(b)./pi)*2*pi;
                b = (st3angles < -pi) | (st3angles > pi);
                st3angles(b) = st3angles(b) - floor(st3angles(b)./pi)*2*pi;
                
                [~,lineNodeIdx(1)] = max(st1angles);
                [~,lineNodeIdx(2)] = min(st1angles);
                [~,lineNodeIdx(3)] = max(st2angles);
                [~,lineNodeIdx(4)] = min(st2angles);
                [~,lineNodeIdx(5)] = max(st3angles);
                [~,lineNodeIdx(6)] = min(st3angles);
                
                A = mod([st1angles(lineNodeIdx(1)) + refangles(1); ...
                    st1angles(lineNodeIdx(2)) + refangles(1); ...
                    st2angles(lineNodeIdx(3)) + refangles(2); ...
                    st2angles(lineNodeIdx(4)) + refangles(2); ...
                    st3angles(lineNodeIdx(5)) + refangles(3); ...
                    st3angles(lineNodeIdx(6)) + refangles(3)], 2*pi);
                
                if A(6) > A(5)
                    check = A(1:4) > A(5) & A(1:4) < A(6);
                else
                    check = (A(1:4) > A(5) & A(1:4) < 2*pi) | (A(1:4) > 0 & A(1:4) < A(6));
                end
                
                if check
                    D = sum(mod([A(6) - A(1), A(2) - A(3), A(4) - A(5)], 2*pi))/3;
                    subtree2rot = D - A(4) + A(5);
                    subtree1rot = D - A(2) + A(3) + subtree2rot;
                    
                    subtree1Pts = [tr.x(subtree1), tr.y(subtree1)]';
                    subtree2Pts = [tr.x(subtree2), tr.y(subtree2)]';
                    refPt = [tr.x(selectedNode); tr.y(selectedNode)];
                    
                    R1 = [cos(subtree1rot) -sin(subtree1rot); sin(subtree1rot) cos(subtree1rot)];
                    R2 = [cos(subtree2rot) -sin(subtree2rot); sin(subtree2rot) cos(subtree2rot)];
                    
                    subtree1NewPts = bsxfun(@plus,R1*(bsxfun(@minus,subtree1Pts,refPt)),refPt);
                    subtree2NewPts = bsxfun(@plus,R2*(bsxfun(@minus,subtree2Pts,refPt)),refPt);
                    
                    tr.x(subtree1) = subtree1NewPts(1,:);
                    tr.y(subtree1) = subtree1NewPts(2,:);
                    tr.x(subtree2) = subtree2NewPts(1,:);
                    tr.y(subtree2) = subtree2NewPts(2,:);
                    
                    tr.ang(subtree1) = tr.ang(subtree1) + subtree1rot;
                    tr.ang(subtree2) = tr.ang(subtree2) + subtree2rot;
                end
            end
        end
        
        tr.ang = mod(tr.ang, 2*pi) + rotation;
        tr.tang = mod(tr.tang, 2*pi) + rotation;
        tmp = [tr.x, tr.y] * [cos(rotation) sin(rotation); -sin(rotation) cos(rotation)];
        tr.x = tmp(:,1);
        tr.y = tmp(:,2);
end

%--------------------------------------------------------------------------
% Rendering graphic objects
fig = figure;
% set figure tag for testing purposes.
set(fig,'Tag','phytreeplot');
h.axes = axes; 

%set axes tag for testing purpose
set(h.axes,'Tag',orientation);
%set render type and rotation in appdata for testing purposes
setappdata(fig,'RenderType',renderType);
setappdata(fig,'Rotation',rotation*180/pi);
hold on;

% setting the axes
if ~circularType
    switch orientation
        case {'left','right'}
            set(h.axes,'YTick',1:numel(tr.terminalNodes),...
                'XDir','normal','Ydir','reverse',...
                'YtickLabel','','YAxisLocation','Right','Fontsize',9)
        case {'bottom','top'}
            set(h.axes,'XTick',1:numel(tr.terminalNodes),...
                'Xdir','normal','Ydir','normal',...
                'XtickLabel','','XAxisLocation','Top','Fontsize',9)
    end
else
    set(h.axes,'XTick',[],'YTick',[],'Fontsize',9)
    set(h.axes,'Position',[.05 .05 .9 .9])
    dispTerminalLabels = false;
    axis equal
end

% drawing lines
switch renderType
    case 'square'
        switch orientation
            case {'left','right'}
                X = tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]);
                Y = tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]);
            case {'top','bottom'}
                Y = tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]);
                X = tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        end
    case 'angular'
        switch orientation
            case {'left','right'}
                X = tr.x([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
                Y = tr.y([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
            case {'top','bottom'}
                Y = tr.x([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
                X = tr.y([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        end
    case 'radial'
        rho = tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]);
        theta = tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        
        offset = pi/40;
        rho = rho - min(rho(:));
        theta = theta - min(theta(:));
        theta = theta*(2*pi - 2*offset)/max(theta(:));
        theta = theta + offset;
        
        n_interp_pts = 50;
        r_interp = zeros(n_interp_pts+1, size(rho,2));
        t_interp = zeros(n_interp_pts+1, size(theta,2));
        r_interp(1,:) = rho(1,:);
        r_interp(2:n_interp_pts+1,:) = [bsxfun(@plus, rho(2,:), bsxfun(@times, (0:n_interp_pts-2)', rho(3,:) - rho(2,:))/(floor(n_interp_pts) - 1)); rho(3,:)];
        t_interp(1,:) = theta(1,:);
        t_interp(2:n_interp_pts+1,:) = [bsxfun(@plus, theta(2,:), bsxfun(@times, (0:n_interp_pts-2)', theta(3,:) - theta(2,:))/(floor(n_interp_pts) - 1)); theta(3,:)];
        
        [X,Y] = pol2cart(t_interp, r_interp);
        
        tmp = [X(:), Y(:)] * [cos(rotation) sin(rotation); -sin(rotation) cos(rotation)];
        X = reshape(tmp(:,1), size(t_interp));
        Y = reshape(tmp(:,2), size(t_interp));
        
        tr.x = X(1,:)';
        tr.y = Y(1,:)';
        tr.ang = theta(1,:)' + rotation;
        tr.tang = theta(1,:)' + rotation;
        
    case {'equalangle', 'equaldaylight'}
        X = [tr.x(tr.par)', 0; tr.x'];
        Y = [tr.y(tr.par)', 0; tr.y'];
end
dataRange = [min(X(:)) max(X(:)) min(Y(:)) max(Y(:))];
h.BranchLines = plot(X,Y,'-k');

% drawing nodes
if ~circularType
    switch orientation
        case {'left','right'}
            h.BranchDots = plot(tr.x(prunedBranchIndex),tr.y(prunedBranchIndex),'o',...
                'MarkerSize',5,'MarkerEdgeColor','k',...
                'MarkerFaceColor','b');
            h.LeafDots = plot(tr.x(prunedLeafIndex),tr.y(prunedLeafIndex),'square',...
                'MarkerSize',4,'MarkerEdgeColor','k',...
                'MarkerFaceColor','w');
        case {'top','bottom'}
            h.BranchDots = plot(tr.y(prunedBranchIndex),tr.x(prunedBranchIndex),'o',...
                'MarkerSize',5,'MarkerEdgeColor','k',...
                'MarkerFaceColor','b');
            h.LeafDots = plot(tr.y(prunedLeafIndex),tr.x(prunedLeafIndex),'square',...
                'MarkerSize',4,'MarkerEdgeColor','k',...
                'MarkerFaceColor','w');
    end
else
    h.BranchDots = plot(tr.x(prunedBranchIndex),tr.y(prunedBranchIndex),'o',...
        'MarkerSize',5,'MarkerEdgeColor','k',...
        'MarkerFaceColor','b');
    h.LeafDots = plot(tr.x(prunedLeafIndex),tr.y(prunedLeafIndex),'square',...
        'MarkerSize',4,'MarkerEdgeColor','k',...
        'MarkerFaceColor','w');
end

% set branch node labels
if ~circularType && any(strcmp(orientation,{'bottom','top'}))
    X = tr.y(prunedBranchIndex);
    Y = tr.x(prunedBranchIndex);
else % {'left', 'right'}
    X = tr.x(prunedBranchIndex);
    Y = tr.y(prunedBranchIndex);
end
h.branchNodeLabels = text(X,Y,tr.names(prunedBranchIndex));
vset(h.branchNodeLabels,'Anchor',[X,Y])
set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
set(h.branchNodeLabels,'interpreter','none')
set(h.branchNodeLabels,'Fontsize',fontSizeBranch);
if ~circularType
    switch orientation
        case {'left','right'}
            set(h.branchNodeLabels,'vertical','bottom')
            set(h.branchNodeLabels,'horizontal','right')
            set(h.branchNodeLabels(~ismember(prunedBranchIndex,branchIndex)),'horizontal','left')
        case {'bottom','top'}
            set(h.branchNodeLabels,'vertical','bottom')
            set(h.branchNodeLabels,'Rotation',30)
    end
else
    hal = {'right','left'};
    val = {'top','bottom'};
    vset(h.branchNodeLabels,'horizontalalignment',hal((X>0)+1)');
    vset(h.branchNodeLabels,'verticalalignment',val((Y>0)+1)');
end


% set leaf nodes labels
if ~circularType && any(strcmp(orientation,{'bottom','top'}))
    X = tr.y(prunedLeafIndex);
    Y = tr.x(prunedLeafIndex);
else % {'left', 'right'}
    X = tr.x(prunedLeafIndex);
    Y = tr.y(prunedLeafIndex);
end
h.leafNodeLabels = text(X,Y,tr.names(prunedLeafIndex));
vset(h.leafNodeLabels,'Anchor',[X,Y])
set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
set(h.leafNodeLabels,'interpreter','none')
set(h.leafNodeLabels,'Fontsize',fontSizeLeaf);
if ~circularType
    switch orientation
        case {'bottom','top'}
            set(h.leafNodeLabels,'Rotation',60)
    end
else
    hal = {'left', 'right'};
    val = {'top','bottom'};
    if rotateLeafLabels && numel(prunedLeafIndex)>0
        rot = tr.ang(prunedLeafIndex)*180/pi;
        
        [~,k] = sort(tr.tang(prunedLeafIndex)*180/pi);
        rot = rot(k);
        rot = cumsum([rem(rot(1)+3600,360);diff(rot)+(diff(rot)<-180).*360-(diff(rot)>180).*360]);
        rot(end+1) = rot(1)+360;
        [d,w] = min(diff(rot));
        while d<-0.1
            rot(w+[0 1]) = mean(rot(w+[0 1])) + range(rot(w+[0 1]))./10.*[-1 1];
            if w==numel(rot)-1
                rot(1) = rot(end)-360;
            elseif w==1
                rot(end) = rot(1)+360;
            end
            [d,w] = min(diff(rot));
        end
        nrot(k) = rot(1:end-1);
        nrot = rem(nrot,360);
        flip = nrot>90 & nrot<270;
        nrot(flip) = nrot(flip) - 180;
        vset(h.leafNodeLabels,'rotation', num2cell(nrot))
        vset(h.leafNodeLabels,'horizontalalignment',hal((flip)+1)');
        vset(h.leafNodeLabels,'verticalalignment',val((flip)+1)');
    else
        vset(h.leafNodeLabels,'horizontalalignment',hal((X<0)+1)');
        vset(h.leafNodeLabels,'verticalalignment',val((Y>0)+1)');
    end
end

% set terminal nodes labels
if ~circularType
    if any(strcmp(orientation,{'bottom','top'}))
        X = tr.y(tr.terminalNodes);
        Y = tr.x(tr.terminalNodes) * 0;
    else % {'left', 'right'}
        X = tr.x(tr.terminalNodes) * 0;
        Y = tr.y(tr.terminalNodes);
    end
    h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes));
    vset(h.terminalNodeLabels,'Anchor',[X,Y])
    switch orientation
        case {'bottom','top'}
            set(h.terminalNodeLabels,'Rotation',90)
    end
else
    %for circular layouts there are not terminal Labels, anyways, we create
    %an empty text to be consistent in the output structure with handles
    h.terminalNodeLabels = text(0,0,' ');
    vset(h.terminalNodeLabels,'Anchor',[0,0])
end
set(h.terminalNodeLabels,'interpreter','none')
set(h.terminalNodeLabels,'Fontsize',fontSizeTerminal);


% turns-off labels accordingly to input parameters
if ~dispTerminalLabels
    set(h.terminalNodeLabels,'visible','off');
end
if ~dispBranchLabels
    set(h.branchNodeLabels,'visible','off');
end
if ~dispLeafLabels
    set(h.leafNodeLabels,'visible','off');
end

box on
hold off

% -------------------------------------------------------------------------
% Adjust tree layout for best rendering on display

if ~circularType
% Minimum distance between the tree extent and the axes edges (set as a
% ratio), this setback is used when the extent of labels do not govern axis
treeExtentToAxesRatio = 1/20;
extDataRange = dataRange+dataRange*[-1 0;1 0;0 -1;0 1]*treeExtentToAxesRatio*[-1 1 0 0;0 0 -1 1];

% Minimum tree width (in Pixels) which is good enough for visualization. 
% When the label extents do not push back, the actual width will be larger.
% Note: Width extent is horizontal in the screen for 'left' and 'right'
% orientation and vertical for 'top' and 'bottom' orientation.
minimumWidthPixels = 100; % (pixels)

% Extra pixels around labels to avoid overlapping of the labels with other
% graph objects
extraPixels = 15;

% Required tree height (in Pixels) to minimize labels overlap.
% Note: Height extent is vertical in the screen for 'left' and 'right'
% orientation and horizontal for 'top' and 'bottom' orientation.
reqHeightPixels = fontSizeBranch * 1.25 * 4/3 * tr.numLeaves + 60;
                      % 1.25 = vertical line spacing
                      % 4/3  = point to pixel conversion
                      % 60   = axes to figure border (in pixels)

set(h.axes,'Units','pixels')
figurePosition = get(fig,'Position');
originalFigurePosition = figurePosition;

% Create a preliminary axes position and axes limits to figure out desired
% label extents and required pixels
if any(strcmp(orientation,{'left','right'}))
    prelimX = minimumWidthPixels;
    prelimY = reqHeightPixels;
else
    prelimX = reqHeightPixels;
    prelimY = minimumWidthPixels;
end
set(h.axes,'Position',[1 1 prelimX prelimY]);
axis(dataRange)
adjustlabelpos([h.branchNodeLabels;h.leafNodeLabels],extraPixels,1)
% find required extents (data units) for all elements inside the axes
if dispLeafLabels 
    if dispBranchLabels
        reqext = getjointextent([h.branchNodeLabels;h.leafNodeLabels],extraPixels);
    else
        reqext = getjointextent(h.leafNodeLabels,extraPixels);
    end
else
    if dispBranchLabels
        reqext = getjointextent(h.branchNodeLabels,extraPixels);
    else
        reqext = extDataRange;
    end
end
reqext = min([reqext.*[1 -1 1 -1];extDataRange.*[1 -1 1 -1]]).*[1 -1 1 -1];
% change required extent to pixels
data2pixel = [prelimX prelimY]./(dataRange*[-1 0;1 0;0 -1;0 1]);
reqextpix = [diff(reqext([1 2])) diff(reqext([3 4]))].*data2pixel;
% find required extents (data units) for all elements inside the figure
reqextpixfig = reqextpix + 60; % border pixels around the axes
if dispTerminalLabels
    reqextpixtl = getjointextent([h.terminalNodeLabels],extraPixels) ...
                  * [-1 0;1 0;0 -1;0 1] .* data2pixel;
    if any(strcmp(orientation,{'left','right'}))
        reqextpixfig(1) =  reqextpixfig(1) + reqextpixtl(1) - 30;
    else % {'bottom','top'}
        reqextpixfig(2) =  reqextpixfig(2) + reqextpixtl(2) - 30;
    end
end
% adjust figure size according the required extent for the HEIGHT of the
% tree and the available space on the screen
if any(strcmp(orientation,{'left','right'}))
    bottomFigPosition = figurePosition(2) + figurePosition(4) - reqextpixfig(2);
    if bottomFigPosition>=1
        if reqextpixfig(2)>figurePosition(4)
            figurePosition(2) = bottomFigPosition;
            figurePosition(4) = reqextpixfig(2);
        end
    else
        figurePosition(2) = 1;
        figurePosition(4) = min(reqextpixfig(2),get(0,'ScreenSize')*[0;0;0;1]-70);
        % 70  = space for the figure toolbar, menubar and titlebar
    end
else % {'bottom','top'}
    rightFigPosition = max(figurePosition(3),reqextpixfig(1))+figurePosition(1);
    rightFigPosition = min(rightFigPosition,get(0,'ScreenSize')*[0;0;1;0]);
    leftFigPosition = max(1,rightFigPosition - reqextpixfig(1));
    leftFigPosition = min(figurePosition(1),leftFigPosition);
    figurePosition(1) = leftFigPosition;
    figurePosition(3) = rightFigPosition-leftFigPosition;
end

% adjust figure size according the required extent for the WIDTH of the
% tree, in the case of the WIDTH we allow the figure limits go beyond the
% limits of the screen
if any(strcmp(orientation,{'left','right'}))    
    if reqextpixfig(1) > figurePosition(3)
        figurePosition(3) = reqextpixfig(1);
    end
else % {'bottom','top'}
    if reqextpixfig(2) > figurePosition(4)
        figurePosition(2) = sum(figurePosition([2,4])) - reqextpixfig(2);
        figurePosition(4) = reqextpixfig(2);
    end
end

% set final figure position and axes position
set(fig,'Position',figurePosition); 

% takes care of special case in which the required figure size is limited
% by the screen but the origin stays in negative coordinates of the screen
obtainedFigPos = get(fig,'Position');
if figurePosition(4)-obtainedFigPos(4)>1
    obtainedFigPos(2) = obtainedFigPos(2) + (figurePosition(4)-obtainedFigPos(4));
    set(fig,'Position',obtainedFigPos); 
end

if any(strcmp(orientation,{'left','right'}))
    if dispTerminalLabels
        minextpix = originalFigurePosition(3)-30-reqextpixtl(1);
    else
        minextpix = originalFigurePosition(3)-60;
    end
    axesPosition = [30 30 max(minextpix,reqextpix(1)) figurePosition(4)-60];
else % {'bottom','top'}
    if dispTerminalLabels
        minextpix = originalFigurePosition(4)-30-reqextpixtl(2);
    else
        minextpix = originalFigurePosition(4)-60;
    end
    axesPosition = [30 30 figurePosition(3)-60 max(minextpix,reqextpix(2))];
end
set(h.axes,'Position',axesPosition)

% This loop readjusts the label positions and then adjust the axis limits
% according the new label extent iteratively until the axis limits are 1%
% at most off the optimal value
reqext = axis;
reqextOld = [inf inf inf inf];
for i =1:10
    if max(abs(reqextOld - reqext)./[diff(xlim) diff(xlim) diff(ylim) diff(ylim)])<.01
        break
    end
    reqextOld = reqext;
    adjustlabelpos([h.leafNodeLabels;h.branchNodeLabels],extraPixels,1)
    % find required extents (data units) for all elements inside the axes
    if dispLeafLabels 
        if dispBranchLabels
            reqext = getjointextent([h.branchNodeLabels;h.leafNodeLabels],extraPixels);
        else
            reqext = getjointextent(h.leafNodeLabels,extraPixels);
        end
    else
        if dispBranchLabels
            reqext = getjointextent(h.branchNodeLabels,extraPixels);
        else
            reqext = extDataRange;
        end
    end
    reqext = min([reqext.*[1 -1 1 -1];extDataRange.*[1 -1 1 -1]]).*[1 -1 1 -1];
    axis(reqext)
end

% once finished moving the axis limits, put the terminal branches in their
% final location along the axes"
tpos = vget(h.terminalNodeLabels,'Anchor');
if any(strcmp(orientation,{'left','right'}))
    tpos(:,1) = max(xlim);
else % {'bottom','top'}
    tpos(:,2) = max(ylim);
end
vset(h.terminalNodeLabels,'Anchor',tpos)
adjustlabelpos(h.terminalNodeLabels,extraPixels,1)

% reverse X axis in case of 'right' orientation and adjust labels
if strcmp(orientation,'right')
    set(h.axes,'Xdir','reverse')
    axesPosition(1) = figurePosition(3)-sum(axesPosition([1 3]));
    set(h.axes,'Position',axesPosition)
    set(h.branchNodeLabels,'horizontal','left')
    set(h.branchNodeLabels(~ismember(prunedBranchIndex,branchIndex)),'horizontal','right')
    set(h.leafNodeLabels,'horizontal','right')
    set(h.terminalNodeLabels,'horizontal','right')
end

% reverse Y axis in case of 'top' orientation and adjust labels
if strcmp(orientation,'top')
    set(h.axes,'Ydir','reverse')
    axesPosition(2) = figurePosition(4)-sum(axesPosition([2 4]));
    set(h.axes,'Position',axesPosition)
    set(h.branchNodeLabels,'horizontal','left')
    set(h.leafNodeLabels,'horizontal','right')
    set(h.terminalNodeLabels,'horizontal','right')
    if dispLeafLabels 
        if dispBranchLabels
            reqext = getjointextentyrev([h.branchNodeLabels;h.leafNodeLabels],extraPixels);
        else
            reqext = getjointextentyrev(h.leafNodeLabels,extraPixels);
        end
    else
        if dispBranchLabels
            reqext = getjointextentyrev(h.branchNodeLabels,extraPixels);
        else
            reqext = extDataRange;
        end
    end
    reqext = min([reqext.*[1 -1 1 -1];extDataRange.*[1 -1 1 -1]]).*[1 -1 1 -1];
    axis(reqext)
    tpos = vget(h.terminalNodeLabels,'Anchor');
    tpos(:,2) = max(ylim);
    vset(h.terminalNodeLabels,'Anchor',tpos)
    adjustlabelpos(h.terminalNodeLabels,extraPixels,-1)
end

set(h.axes,'Units','Normalized')

else % circularType

set(h.axes,'Units','Pixels')
figurePosition = get(fig,'Position');
originalFigurePosition = figurePosition;

% Readjust dataRange to reflect an equal axes
dataRange = dataRange*[.5 .5 0 0;.5 .5 0 0;0 0 .5 .5;0 0 .5 .5] + ...
            max(dataRange*[-.5 0;.5 0;0 -.5;0 .5])*[-1 1 -1 1];

% Minimum distance between the tree extent and the axes edges (set as a
% ratio), this setback is used when the extent of labels do not govern axis
treeExtentToAxesRatio = 1/10;
extDataRange = dataRange+dataRange*[-1 0;1 0;0 -1;0 1]*treeExtentToAxesRatio*[-1 1 0 0;0 0 -1 1];

% Minimum tree width (in Pixels) which is good enough for visualization. 
% When the label extents do not push back, the actual width will be larger.
% Note: Width extent is horizontal in the screen for 'left' and 'right'
% orientation and vertical for 'top' and 'bottom' orientation.
minimumWidthPixels = 250; % (pixels)    

% Create a preliminary axes position and axes limits to figure out desired
% label extents and required pixels
set(h.axes,'Position',[1 1 minimumWidthPixels minimumWidthPixels]);
axis(dataRange)

% find required extents (data units) for all elements inside the axes
if dispLeafLabels 
    if dispBranchLabels
        reqext = getjointextent([h.branchNodeLabels;h.leafNodeLabels],0);
    else
        reqext = getjointextent(h.leafNodeLabels,0);
    end
else
    if dispBranchLabels
        reqext = getjointextent(h.branchNodeLabels,0);
    else
        reqext = extDataRange;
    end
end
reqext = min([reqext.*[1 -1 1 -1];extDataRange.*[1 -1 1 -1]]).*[1 -1 1 -1];

% change required extent to pixels
data2pixel = minimumWidthPixels./diff(dataRange([1 2]));
reqextpix = [diff(reqext([1 2])) diff(reqext([3 4]))].*data2pixel;

maxPixelsAvailable = get(0,'ScreenSize')*[0 0;0 0;1 0;0 1]-[60 130];

finextpix = max([min([maxPixelsAvailable;reqextpix]);originalFigurePosition([3 4])-60]);

figurePosition([3 4]) = [finextpix(1)+60 finextpix(2)+60];
figurePosition(2) = sum(originalFigurePosition([2 4]))-finextpix(2)-60;

if figurePosition(1)+figurePosition(3)>maxPixelsAvailable(1)
    figurePosition(1) = maxPixelsAvailable(1)-figurePosition(3)+60;
end

figurePosition(2) = max(figurePosition(2),1);

set(fig,'Position',figurePosition); 

set(h.axes,'Position',[30 30 finextpix(1) finextpix(2)])

for i = 1:15
    if dispLeafLabels
        if dispBranchLabels
            reqext = getjointextent([h.branchNodeLabels;h.leafNodeLabels],0);
        else
            reqext = getjointextent(h.leafNodeLabels,0);
        end
    else
        if dispBranchLabels
            reqext = getjointextent(h.branchNodeLabels,0);
        else
            reqext = extDataRange;
        end
    end
    reqext = min([reqext.*[1 -1 1 -1];extDataRange.*[1 -1 1 -1]]).*[1 -1 1 -1];
    if any(reqext*[-1 0;1 0;0 -1;0 1]>finextpix./data2pixel)
        axis(reqext)
        break
    end
    axis(reqext)
end
set(h.axes,'Units','Normalized')
end
    
% -------------------------------------------------------------------------    
% store handles in 'UserData', this is not the preferred because now the
% handles are returned in the first input argument, however we keep it for
% backwards compatibility.
set(fig,'UserData',h)
if nargout
    handles = h;
end

% we do not return the handles unless the user specifically ask for them
if nargout==0
    clear h
end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function extraXYdis = getExtraXYdis(extraPixels)
% Calculates extra distance that needs to be added to the labels to avoid
% overlaying on other graph objects.
% Note: assumes the current units of axes are 'pixels'
% if ~strcmpi('Pixels',get(gca,'Units'))
%     error('Bioinfo:phytree:plot:getExtraXYdis:InvalidAxesUnits',...
%           'Invalid Axes Units');
% end
extraXYdis = extraPixels./get(gca,'Pos')*[0 0;0 0;diff(xlim) 0;0 diff(ylim)];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adjustlabelpos(h,extraPixels,ydir)
XY = getExtraXYdis(extraPixels)./2;
pos = vget(h,'Anchor');
theta = vget(h,'Rotation')*pi./180;
theta = theta + strcmp(vget(h,'HorizontalAlignment'),'right')*pi;
pos = [pos(:,1)+XY(1).*cos(theta) pos(:,2)+ydir.*XY(2).*sin(theta) theta.*0];
vset(h,'Position',pos)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ext = getjointextent(h,extraPixels)
XY = getExtraXYdis(extraPixels);
theta = vget(h,'Rotation')*pi./180;
ext = vget(h,'Extent');
ext = [ext(:,1)-abs(XY(1).*cos(theta))./2 ext(:,2)-abs(XY(2).*sin(theta))./2 ...
       ext(:,3)+abs(XY(1).*cos(theta))    ext(:,4)+abs(XY(2).*sin(theta))];
ext = [min(ext(:,1)) max(sum(ext(:,[1,3]),2)) ...
       min(ext(:,2)) max(sum(ext(:,[2,4]),2))];
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ext = getjointextentyrev(h,extraPixels)
XY = getExtraXYdis(extraPixels);
theta = vget(h,'Rotation')*pi./180;
ext = vget(h,'Extent');
ext = [ext(:,1)-abs(XY(1).*cos(theta))./2 ext(:,2)+abs(XY(2).*sin(theta))./2 ...
       ext(:,3)+abs(XY(1).*cos(theta))    ext(:,4)+abs(XY(2).*sin(theta))];
ext = [min(ext(:,1)) max(sum(ext(:,[1,3]),2)) ...
       min(ext(:,2)-ext(:,4)) max(ext(:,2)) ];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = vget(h,property)
% Get a property or application_data of a group handles and put it into a
% numeric matrix when possible, handles the special case of one handle.
switch numel(h)
    case 0
        val = {};
    case 1
        if isprop(h,property)
            val = {get(h,property)};
        else
            val = {getappdata(h,property)};
        end
    otherwise
        if all(isprop(h,property))
            val = get(h,property);
        else
            val = cell(numel(h),1);
            for i = 1:numel(h)
                val{i} = getappdata(h(i),property);
            end
        end
end
if all(cellfun(@(x) isvector(x)&&isnumeric(x),val))
    val = cell2mat(val);
end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vset(h,property,val)
% Set a property or application_data of a group handles 
for i = 1:numel(h)
    if iscell(val)
        v = val{i};
    else
        v = val(i,:);
    end
    if isprop(h(i),property)
        set(h(i),property,v)
    else
        setappdata(h(i),property,v)
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = getChildNodes(node, tr, varargin)

if isempty(varargin)
    in = [];
elseif numel(varargin) == 1
    in = varargin{1};
end

out = [in node];
    
if node > tr.numLeaves
    out = getChildNodes(tr.tree(node - tr.numLeaves, 1), tr, out);
    out = getChildNodes(tr.tree(node - tr.numLeaves, 2), tr, out);
end
end
