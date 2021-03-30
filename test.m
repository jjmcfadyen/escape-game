maze = zeros(MDP.dim(1)*2-1,MDP.dim(2)*2-1);

% add inner borders
for r = 2:2:size(maze,1)
    maze(r,:) = NaN; 
end
for c = 2:2:size(maze,2)
    maze(:,c) = NaN; 
end

% ignore diagonal
maze(2:2:end,2:2:end) = -Inf;

% randomly connect the rooms until all are connected
M = maze;
generating = true;
while generating
   
    % open a random connection
    nan_states = find(isnan(M));
    disp(['Wall count = ' num2str(length(nan_states))])
    M(nan_states(randi(length(nan_states)))) = 0;
    
    % check if all states are connected
    states = find(M==0);
    [I,J] = find(M==0);
    coords = [I,J];
    
    A = pdist2(coords,coords) == 1;

    for st = 1:length(states)
        tree = states(st);
        searching = true;
        while searching
            cc = length(tree);
            for i = 1:length(tree)
                next = setdiff(states(find(A(states==tree(i),:))),tree);
                tree = [tree; setdiff(next,tree)];
            end
            if length(tree) == cc
                searching = false;
            end
        end
        if all(ismember(states,unique(tree)))
            generating = false;
            break
        end
    end
    
end

% add outside border walls
maze(1,:) = NaN;
maze(end,:) = NaN;
maze(:,1) = NaN;
maze(:,end) = NaN;

