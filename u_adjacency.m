function MDP = u_adjacency(MDP)
% --------------------------------------------------------
% Converts a map of state numbers into an adjacency matrix
% --------------------------------------------------------
%% Transition matrix (nStates x nStates) - from x to

% Identify coordinates of viable states
[I,J] = find(MDP.map);
rr = [J,I]; % row, col indices of each state

MDP.viableStates = sort(MDP.map(~isnan(MDP.map)));
MDP.nViable = length(MDP.viableStates);

% Compute adjacency matrix
A = pdist2(rr,rr) == 1;

% Define one-step transition matrix
T = double(A);
T(T==0) = NaN;

% Update terminal states
if isfield(MDP,'terminal')
    if ~isempty(MDP.terminal)
        for t = 1:length(MDP.terminal)
            T(MDP.terminal(t),:) = NaN;
            T(MDP.terminal(t),MDP.terminal(t)) = 1;
        end
    end
end

% Allow transitions to same state (i.e. self-transitions)
T(find(eye(size(T,1),size(T,2)))) = 1;

% Insert NaNs
nanidx = find(isnan(MDP.map')); % transpose so that state numbers read left to right, starting from the top
T(nanidx,:) = NaN;
T(:,nanidx) = NaN;

% Insert into MDP variable
MDP.T = T;

end