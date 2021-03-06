function [path,randomness] = mdp_astar(A,B,MDP)

% get start and end coordinates
if length(A)==1
    aidx = u_coordswitch(A,MDP);
else
    aidx = A;
end
if length(B)==1
    bidx = u_coordswitch(B,MDP);
else
    bidx = B;
end

% create binary map (excluding starting state)
% --- 1 = walkable
% --- 0 = obstacle
map = MDP.map;
map(~isnan(map)) = 1;
map(isnan(map)) = 0;

map(aidx(1),aidx(2)) = 0;

% create empty g/f/h cost maps
tmp = map;
tmp(tmp==0) = NaN;
tmp(tmp==1) = Inf;

gcost = tmp;
fcost = tmp;
hcost = tmp;

% create map of which states are open/closed
% --- 1 = open
% --- 0 = closed
omap = zeros(size(map,1),size(map,2));

% create map of parent nodes
parents = nan(size(map,1),size(map,2));

% iteratively open up the map
position = aidx;
closed = nan(0,2);
randomness = false;
alternatives = []; % position, options left
while true
    
    % open up the surrounding states
    coords = [];
    for a = 1:4 % north, east, south, west
        coord = position;
        if a == 1
            coord(1) = coord(1)-1;
        elseif a == 2
            coord(2) = coord(2)+1;
        elseif a == 3
            coord(1) = coord(1)+1;
        elseif a == 4
            coord(2) = coord(2)-1;
        end
        if all(coord>0) && coord(1) < size(map,1) && coord(2) < size(map,2)
            if map(coord(1),coord(2)) == 1
                omap(coord(1),coord(2)) = 1;
                coords = [coords; coord];
            end
        end
    end
    
    % calculate their g/h/f costs
    g = pdist2(aidx,coords);
    h = pdist2(bidx,coords);
    f = g + h;
    
    newg = gcost;
    newh = hcost;
    newf = fcost;
    
    for c = 1:size(coords,1)
        I = coords(c,1);
        J = coords(c,2);
        if isempty(closed) || ~any(ismember(closed,[I J],'rows'))
            newg(I,J) = g(c);
            newh(I,J) = h(c);
            newf(I,J) = f(c);
        end
    end
    
    % only update if the f values are an improvement
    idx = (fcost - newf) > 0;
    fcost(idx) = newf(idx);
    gcost(idx) = newg(idx);
    hcost(idx) = newh(idx);
    
    parents(idx) = u_coordswitch(position,MDP);
    
%{     
    pmap = fcost;
    pmap(pmap==Inf) = 0;

    imagesc(pmap); hold on
    cmap = [0 0 0; u_colours(max(pmap(:))-min(pmap(:))+1,'viridis')];
    colormap(cmap);
    caxis([min(pmap(:))-1 max(pmap(:))])
    
    scatter(position(end,2),position(end,1),100,'markerfacecolor',[.5 .5 .5],'markeredgecolor','k'); hold on
    
    for c = 1:size(closed,1)
        scatter(closed(c,2),closed(c,1),100,'x','r'); hold on
    end
%}
    if ismember(position,bidx,'rows')
        
        % trace back the optimal path
        current = position;
        path = MDP.map(position(1),position(2));
        while true
            idx = parents(current(1),current(2));
            if A == idx
                break
            end
            path = [idx; path];
            current = u_coordswitch(path(1),MDP);
        end
        
        break
    end
    
    % decide which state to close, based on entier map of g/f/h costs
    [I,J] = find(fcost ==  min(fcost(:)));
    min_f = [I,J];
    
    [I,J] = find(hcost ==  min(hcost(:)));
    min_h = [I,J];
    
    % if there is only ONE node to close...
    if size(min_f,1) == 1
        cnode = min_f; 
    elseif size(min_h,1) == 1
        cnode = min_h;
    else
        rand_choice = randi(size(min_f,1));
        cnode = min_f(rand_choice,:);
        warning('There might be two shortest paths')
        randomness = true;
        unchosen = setdiff(1:size(min_f,1),rand_choice);
        choices_left = nan(1,4);
        choices_left(1:length(unchosen)) = unchosen;
        alternatives = [alternatives; position choices_left];
    end
    
    % close it
    omap(cnode(1),cnode(2)) = 0;

    % update g/h/f costs
    gcost(cnode(1),cnode(2)) = Inf;
    hcost(cnode(1),cnode(2)) = Inf;
    fcost = gcost + hcost;

    % update closed & parents lists
    closed = [closed; cnode];
    
    % update position
    position = cnode;
    
end

% Try again with alternatives
disp('')

end