function MDP = u_maze(MDP,type)

%% Parameters (automatic generation)

% Obstacles
nSolutions = 2; % number of possible routes from start to finish
odensity   = 0.75; % 1 = keep all obstacles, 0 = remove all obstaces, 0-1 remove portion of obstacles

% Predator spawn point
depth   = 3; % how many moves away (max) from the shortest path to generate the predator (Inf = pick any state in the map)
barrier = 2;   % how many moves away (min) from the agent

%% Set up

if strcmp(type,'interactive')
    generateMap = false; % just randomise the predator location
elseif strcmp(type,'auto')
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
    
    %% Assign safety spots
   
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
    B = u_coordswitch(MDP.terminal,MDP);
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
    
else
    
    %% Manual map definition
    
    % Get larger maze
    maze = zeros(MDP.dim(1)*2+1,MDP.dim(2)*2+1);

    % add outside border walls
    maze(1,:) = NaN;
    maze(end,:) = NaN;
    maze(:,1) = NaN;
    maze(:,end) = NaN;
    
    % Define features
    figure
    for c = 1:5 % walls, start, predator, reward(s), safety
        
        if c < 5 || (c == 5 && MDP.safety.on)
            clicking = true;
            while clicking

                clf
                imagesc(maze); hold on
                for i = 2:size(maze,1)-1
                    for j = 2:size(maze,1)-1
                        scatter(i,j,50,'marker','s','markerfacecolor','k','markeredgecolor','k'); hold on
                    end
                end

                if c == 1
                    scatter(1,size(maze,1),100,'markerfacecolor','g','markeredgealpha',0); hold on
                    scatter(size(maze,2),size(maze,1),100,'markerfacecolor','r','markeredgealpha',0); hold on
                    disp('::: Click on squares to convert to walls and press enter')
                    disp('::: Click green and enter when finished, or red and enter to start again')
                end

                if c > 2
                    scatter(MDP.start(2),MDP.start(1),75,'markerfacecolor','b','markeredgecolor','k'); hold on
                end
                if c > 3
                    scatter(MDP.punishment(2),MDP.punishment(1),75,'markerfacecolor','r','markeredgecolor','k'); hold on
                end
                if c > 4
                    for i = 1:size(MDP.reward,1)
                         scatter(MDP.reward(i,2),MDP.reward(i,1),75,'markerfacecolor',[0 1 0],'markeredgecolor','k'); hold on
                    end
                end

                if c > 1
                    if c == 2
                        title('Choose start location')
                    elseif c == 3
                        title('Choose predator location')
                    elseif c == 4
                        title('Choose reward locations')
                    elseif c == 5
                        title('Choose safe locations')
                    end
                    disp('::: Click on squares to convert and press enter')
                end

                colormap('gray')

                [x,y] = ginput;
                x = round(x);
                y = round(y);

                if c == 1
                    
                    for i = 1:length(x)
                        if isnan(maze(y(i),x(i)))
                            maze(y(i),x(i)) = 0;
                        elseif maze(y(i),x(i)) == 0
                            maze(y(i),x(i)) = NaN;
                        end
                    end

                    if any(ismember([x y],[1 size(maze,1)],'rows'))
                        maze(end,1) = NaN;
                        clicking = false; % finished
                    end

                    if any(ismember([x y],[size(maze,1) size(maze,1)],'rows'))
                        disp('RESETTING')
                        maze = zeros(size(maze,1),size(maze,2));
                        maze(1,:) = NaN;
                        maze(end,:) = NaN;
                        maze(:,1) = NaN;
                        maze(:,end) = NaN;
                    end
                    
                else
                    if c < 4
                        if length(x) > 1
                            warning('Too many locations chosen. Only choose ONE.')
                        elseif isnan(maze(y,x))
                            warning('Must choose a state that is not a wall')
                        else
                            if c == 2
                                MDP.start = [y x];
                            elseif c == 3
                                MDP.punishment = [y x];
                            end
                            clicking = false;
                        end
                    elseif c == 4
                        MDP.reward = [];
                        for i = 1:length(x)
                            if isnan(maze(y(i),x(i)))
                                warning('Reward location must not be a wall. Ignoring')
                            else
                                MDP.reward = [MDP.reward; y(i) x(i)];
                            end
                        end
                        clicking = false;
                    elseif c == 5
                        MDP.safety.loc_safe = [];
                        for i = 1:length(x)
                            if isnan(maze(y(i),x(i)))
                                warning('Safe location must not be a wall. Ignoring')
                            else
                                MDP.safety.loc_safe = [MDP.safety.loc_safe; y(i) x(i)];
                            end
                        end
                        clicking = false;
                        for i = 1:size(MDP.safety.loc_safe,1)
                            scatter(MDP.safety.loc_safe(i,2),MDP.safety.loc_safe(i,1),75,'markerfacecolor',[0 1 1],'markeredgecolor','k'); hold on
                        end
                    end
                end
            end
        end
    end
    
    % number the non-nan states
    MDP.map = reshape(1:numel(maze),size(maze,1),size(maze,2))';;
    MDP.map(isnan(maze)) = NaN;
    
    % update dimensions
    MDP.dim = size(MDP.map);
    MDP.nStates = numel(MDP.map);
    
    % update feature locations
    MDP.start = u_coordswitch(MDP.start,MDP);
    MDP.punishment = u_coordswitch(MDP.punishment,MDP);
    
    for i = 1:size(MDP.reward,1)
        MDP.reward(i,1) = u_coordswitch(MDP.reward(i,:),MDP);
    end
    MDP.reward = MDP.reward(:,1);
    
    MDP.terminal = [MDP.punishment MDP.reward(end)]; % only the last reward is the terminal
    
    if MDP.safety.on
        for i = 1:size(MDP.safety.loc_safe,1)
            MDP.safety.loc_safe(i,1) = u_coordswitch(MDP.safety.loc_safe(i,:),MDP);
        end
        MDP.safety.loc_safe = MDP.safety.loc_safe(:,1);
    end
    
    % recompute adjacency matrix
    MDP = u_adjacency(MDP);
    
end

end