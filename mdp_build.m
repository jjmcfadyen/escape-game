function MDP = mdp_build(MDP)

%% Transition matrix

MDP = u_adjacency(MDP);

%% Reward matrix (nStates x nStates)

MDP.R = zeros(MDP.nStates,MDP.nStates);

if MDP.safety.on && ~isempty(MDP.safety.found)
    switch MDP.safety.method
        case 'merged'
            for i = 1:length(MDP.safety.found)
                MDP.R(:,MDP.safety.found(i)) = MDP.safety.sval;
%                 MDP.R(MDP.safety.found(i),MDP.safety.found(i)) = MDP.safety.sval;
            end
    end
end

if isfield(MDP,'reward')
    MDP.R(:,MDP.reward(1)) = MDP.reward(2);
end
if isfield(MDP,'punishment')
    MDP.R(:,MDP.punishment(1)) = MDP.punishment(2);
end

for t = 1:length(MDP.terminal)
    MDP.R(MDP.terminal(t),MDP.terminal(t)) = 0;
end

A = zeros(MDP.nStates,MDP.nStates);
A(MDP.T==1) = 1;

MDP.R = MDP.R .* A;
MDP.R(isnan(MDP.T)) = NaN;

%% Log info of each s' (value & transition probability)

% Find maximum no. of transitions from any state
max_transitions = 5; % north, south, east, west, stay/self (in no particular order)

data = array2table(cell(MDP.nStates,4),'variablenames',{'st','sprime','T','R'});
for st = 1:MDP.nStates
    
    % how many states does this state transition to (including itself)?
    ntrans = sum(~isnan(MDP.T(st,:)));
    
    sprime = nan(1,max_transitions);
    T = nan(1,max_transitions);
    R = nan(1,max_transitions);
    if ntrans > 0
        
        % states that can be transitioned to
        sprime(1:ntrans) = find(~isnan(MDP.T(st,:)));

        % probability of each transition (normalised)
        T(1:ntrans) = 1-MDP.epsilon;
        T(ntrans+1:end) = MDP.epsilon;

        % reward outcome at each transition
        R(1:ntrans) = MDP.R(st,sprime(1:ntrans));

    end
    
    data.st{st} = st;
    data.sprime{st} = sprime;
    data.T{st} = T;
    data.R{st} = R;
    
end

MDP.info = data;

% vN = {'S','Sprime','R','T'};
% MDP.info = array2table(cell(0,length(vN)),'variablenames',vN);
% for sv = 1:MDP.nViable
%     
%     st = MDP.viableStates(sv);
%     
%     st_prime = find(~isnan(MDP.T(st,:)));
%     r = MDP.R(st,st_prime);
%     
%     t = [1-MDP.epsilon, ones(1,length(r)-1)*MDP.epsilon];
%     t = t / sum(t); % normalise
%     
%     arr = array2table(cell(length(st_prime),length(vN)),'variablenames',vN);
%     for sp = 1:length(st_prime)
%         arr.S(sp) = {st};
%         arr.Sprime(sp) = {circshift(st_prime,sp-1)};
%         arr.R(sp) = {circshift(r,sp-1)};
%         arr.T(sp) = {t};
%     end
%     
%     MDP.info = [MDP.info; arr];
%     
% end
% 
% MDP.info.S = cell2mat(MDP.info.S);
% 
% % Make terminal states deterministic
% for t = 1:length(MDP.terminal)
%     MDP.info.T(MDP.info.S == MDP.terminal(t)) = {1}; 
% end

end