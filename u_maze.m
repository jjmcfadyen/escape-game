function MDP = u_maze(MDP,type)

%% Parameters

% Obstacles
nSolutions = 3; % number of possible routes from start to finish
odensity   = 0.75; % 1 = keep all obstacles, 0 = remove all obstaces, 0-1 remove portion of obstacles

% Predator spawn point
depth   = 3; % how many moves away (max) from the shortest path to generate the predator (Inf = pick any state in the map)
barrier = 2;   % how many moves away (min) from the agent

%% Set up

if strcmp(type,'predator')
    generateMap = false; % just randomise the predator location
elseif strcmp(type,'all')
    generateMap = true;  % generate entirely new map + predator
end

if generateMap
    
    % Get adjacency matrix
    MDP = u_adjacency(MDP);
    
    % Pick a random start location along bottom half
    MDP.start = MDP.map(floor(MDP.dim(1)/2)+1:end,:);
    MDP.start = randsample(MDP.start(:),1);

    % Pick a random end location along top row
    MDP.terminal = MDP.map(1:floor(MDP.dim(1)/2),:);
    MDP.terminal = setdiff(MDP.terminal(:),find(MDP.T(MDP.start,:)==1));
    MDP.terminal = randsample(MDP.terminal(:),1);

    % Re-compute adjacency matrix
    MDP = u_adjacency(MDP);

    %% Randomly walk through

    transitions = [];

    % Get one successful walk from start to end
    completed = false;
    cc = 0;
    while ~completed

        cc = cc + 1;

%         H = figure;
%         bg = checkerboard(1,MDP.dim(1),MDP.dim(2));
%         bg = bg(1:MDP.dim(1),1:MDP.dim(2));
% 
%         imagesc(bg); hold on
%         colormap([1 1 1; .9 .9 .9])
%         set(gca,'ydir','reverse')
% 
%         A = u_coordswitch(MDP.start,MDP);
%         B = u_coordswitch(MDP.terminal,MDP);
%         scatter(A(2),A(1),'markerfacecolor','b','markeredgealpha',0); hold on
%         scatter(B(2),B(1),'markerfacecolor','g','markeredgealpha',0); hold on
% 
%         title(['Iteration ' num2str(cc)])
%         xlim([0.5 MDP.dim(1)+.5])
%         ylim([0.5 MDP.dim(2)+.5])

        visited = [];
        position = MDP.start;

        while true

            visited = [visited; position];

            if any(visited == MDP.terminal)
                completed = true;
                break
            end

            nodes = setdiff(find(MDP.T(position,:)==1),visited);
            if isempty(nodes)
                close(gcf)
                break
            end
            nextnode = nodes(randi(length(nodes)));

%             A = u_coordswitch(position,MDP);
%             B = u_coordswitch(nextnode,MDP);
%             plot([A(2) B(2)],[A(1) B(1)],'b','linewidth',2); hold on
%             drawnow

            position = nextnode;

        end
    end
    routes = {visited};
    transitions = [transitions; visited(1:end-1) visited(2:end)];
    
%     close(H)

    % Take successful path and fork off alternate solutions (from first 2/3 of path, so that it's away from the end)
    if nSolutions > 1
        for s = 1:nSolutions-1

            completed = false;
            cc = 0;
            while ~completed

                cc = cc + 1;

                visited = [];
                position = routes{1}(1:round(length(routes{1})*(2/3)));
                position = position(randi(length(position)));

                while true

                    visited = [visited; position];

                    nodes = setdiff(find(MDP.T(position,:)==1),visited);
                    if isempty(nodes)
                        break
                    end

                    % pick next node, preferencing nodes that aren't in other paths
                    if length(nodes) > 1
                        W = double(~ismember(nodes,unique(transitions(:))));
                        W(W==1) = .75;
                        W(W==0) = .25;
                        nextnode = randsample(nodes,1,true,W);
                    else
                        nextnode = nodes;
                    end

                    position = nextnode;

                    if any(position == MDP.terminal)
                        visited = [visited; position];
                        if ~all(ismember(unique(transitions,'rows'),[visited(1:end-1) visited(2:end)],'rows'))
                            completed = true;
                            break
                        else
                           warning('check') 
                        end
                    end

                end
            end

            % Save successful route
            transitions = [transitions; visited(1:end-1) visited(2:end)];
            routes{s+1} = visited;

%             % Draw
%             for t = 1:length(visited)-1
%                 A = u_coordswitch(visited(t),MDP);
%                 B = u_coordswitch(visited(t+1),MDP);
%                 plot([A(2) B(2)],[A(1) B(1)],'b','linewidth',2); hold on
%                 drawnow
%             end

        end
    end

    % Fill in missing areas
    used = unique(cell2mat(routes'));
    available = setdiff(1:MDP.nStates,used);

    while ~isempty(available)

        position = available(randi(length(available)));
        visited = [];
        while true

            visited = [visited; position];

            if any(visited(end) == used)
                available = setdiff(available,visited);
                used = [used; visited];
                transitions = [transitions; visited(1:end-1) visited(2:end)];
                break
            end

            neighbours = setdiff(find(MDP.T(position,:)==1),visited);
            if isempty(neighbours)
                break;
            end

            best = neighbours(ismember(neighbours,available)); % prioritise going to an empty space
            if length(best) > 0
                nextnode = best(randi(length(best)));
            else
                nextnode = neighbours(randi(length(neighbours)));
            end

%             A = u_coordswitch(position,MDP);
%             B = u_coordswitch(nextnode,MDP);
%             plot([A(2) B(2)],[A(1) B(1)],'b','linewidth',2); hold on
%             drawnow

            position = nextnode;

        end

    end

    transitions = unique(transitions,'rows');

    %% Convert to array

    maze = zeros(MDP.dim(1)*2+1,MDP.dim(2)*2+1);

    % add outside border walls
    maze(1,:) = NaN;
    maze(end,:) = NaN;
    maze(:,1) = NaN;
    maze(:,end) = NaN;

    % add inner borders
    for r = 3:2:size(maze,1)
        maze(r,:) = NaN; 
    end
    for c = 3:2:size(maze,2)
        maze(:,c) = NaN; 
    end

    % add state names
    cc = 0;
    for r = 2:2:size(maze,1)
        cc = cc + 1;
        maze(r,2:2:end) = MDP.map(cc,:);
    end

    % add connections
    for t = 1:size(transitions,1)

        [I,J] = find(maze == transitions(t,1));
        A = [I,J];

        [I,J] = find(maze == transitions(t,2));
        B = [I,J];

        D = (A(:) + B(:)).'/2;

        maze(D(1),D(2)) = 0;

    end

    % clear out obstacles
    if odensity < 1
        innermaze = maze(2:end-1,2:end-1);
        obstacles = find(isnan(innermaze));
        obstacles = randsample(obstacles, round((1-odensity)*length(obstacles)));
        innermaze(obstacles) = 0;
        maze(2:end-1,2:end-1) = innermaze;
    end
    
    MDP.map = reshape(1:numel(maze),size(maze,1),size(maze,2))';
    MDP.map(isnan(maze)) = NaN;
    MDP.dim = size(MDP.map);
    MDP.nStates = numel(maze);

    MDP.start = MDP.map(find(maze==MDP.start));
    MDP.terminal = MDP.map(find(maze==MDP.terminal));

    MDP = u_adjacency(MDP);

end

%% Assign safety spots

if MDP.safety.sdensity > 0
   
    % determine number of safe places
    viable = MDP.map(~isnan(MDP.map));
    nViable = length(viable);
    nSafe = round(MDP.safety.sdensity*nViable);
    
    % evenly distribute across map
    I = u_grid(nSafe);
    
    rows = round(linspace(1,MDP.dim(1),I(1)+1));
    cols = round(linspace(1,MDP.dim(2),I(2)+1));
    
    loc_safe = [];
    for r = 1:length(rows)-1
        for c = 1:length(cols)-1
            innermaze = MDP.map(rows(r):rows(r+1),cols(c):cols(c+1));
            innermaze(ismember(innermaze,MDP.terminal)) = NaN;
            viable = innermaze(~isnan(innermaze));
            if ~isempty(loc_safe)
                viable = setdiff(viable,find(sum(MDP.T(:,loc_safe)==1,2)>0));
            end
            loc_safe = [loc_safe, viable(randi(length(viable)))];
        end
    end
    
    MDP.safety.loc_safe = loc_safe;
    
end

%% Assign predator to random location near shortest path

if isfield(MDP,'reward') || isfield(MDP,'punishment')
    terminal = MDP.reward(1); 
else
    terminal = MDP.terminal; 
end

% get shortest path using A* algorithm
sp = mdp_astar(MDP.start,terminal,MDP);
sp = [MDP.start; sp];

% for p = 1:length(sp)-1
%     [I,J] = find(MDP.map == sp(p));
%     [I(2),J(2)] = find(MDP.map == sp(p+1));
%     plot(J,I,'w','linewidth',2); hold on
% end

% Get random predator location, near to shortest path
if depth == Inf
    nearnodes = setdiff(MDP.map(~isnan(MDP.map)),sp);
else
    for d = 1:depth
        if d == 1
            nearnodes = [];
            for p = 1:length(sp)
                nearnodes = [nearnodes; find(MDP.T(sp(p),:)==1)']; 
            end
        else
            for p = 1:length(nearnodes)
                nearnodes = [nearnodes; find(MDP.T(nearnodes(p),:)==1)']; 
            end
        end
        nearnodes = setdiff(unique(nearnodes(:)),sp);
    end
end

for d = 1:barrier
    if d == 1
        agentbarrier = find(MDP.T(MDP.start,:)==1)';
    else
        for p = 1:length(agentbarrier)
            agentbarrier = [agentbarrier; find(MDP.T(agentbarrier(p),:)==1)'];
        end
    end
end
agentbarrier = unique(agentbarrier);
nearnodes = setdiff(nearnodes,agentbarrier);

MDP.punishment(1) = nearnodes(randi(length(nearnodes)));

% plot
binmaze = ~isnan(MDP.map);

figure;
imagesc(binmaze); hold on
colormap('gray')

set(gca,'ydir','reverse')

A = u_coordswitch(MDP.start,MDP);
B = u_coordswitch(terminal,MDP);
scatter(A(2),A(1),80,'markerfacecolor','b','markeredgecolor','k'); hold on
scatter(B(2),B(1),80,'markerfacecolor','g','markeredgecolor','k'); hold on

[I,J] = find(MDP.map == MDP.punishment(1)); 
scatter(J,I,80,'markerfacecolor','r','markeredgecolor','k')

if isfield(MDP,'safety')
    for s = 1:length(MDP.safety.loc_safe)
        [I,J] = find(MDP.map == MDP.safety.loc_safe(s)); 
        scatter(J,I,80,'markerfacecolor','y','markeredgecolor','k'); hold on
    end
end

end