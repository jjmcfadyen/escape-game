function MDP = mdp_valueIteration(MDP,makeplots)
%     Q-value iteration algorithm.
%
%     
%     Parameters
%     ----------
%     policy : max | min | softmax | pessimism (default = pessimism)
%         Learning rule.
%     gamma : float (default = 0.9)
%         Temporal discounting factor.
%     beta : float (default = 10.0)
%         Inverse temperature for future choice (ignored if policy not softmax).
%     w : float (default = 1.0)
%         Pessimism weight (ignored if policy not pessimism).
%     tol : float, default: 1e-4
%         Tolerance for stopping criteria.
%     max_iter : int, default: 100
%         Maximum number of iterations taken for the solvers to converge.
%
%     References
%     ----------
%     1. Sutton, R. S., & Barto, A. G. (2018). Reinforcement learning: An introduction. MIT press.
%     

%% Defaults

policy = MDP.policy;
gamma = MDP.gamma;
beta = MDP.beta;
w = MDP.w;
tol = MDP.tol;
max_iter = MDP.max_iter;

if isempty(MDP)
    error('MDP needs to be provided')
end
if isempty(policy)
    policy = 'pessimism';
elseif ~any(strcmp(policy,{'max','min','softmax','pessimism','mean'}))
    error(['That policy ("' policy '") does not exist']);
end
if isempty(gamma)
    gamma = 0.95;
elseif gamma < 0 || gamma > 1
    error('Gamma must be between 0 and 1');
end
if isempty(beta)
    beta = 10; % should be less than 50
end
if isempty(w)
    w = 1;
elseif w < 0 || w > 1
    error('w must be between 0 and 1');
end
if isempty(tol)
    tol = 0.0001;
end
if isempty(max_iter)
    max_iter = 100;
end

%% Solve Q values

% Initialise Q values
Q            = cell2mat(MDP.info.sprime);
Q(~isnan(Q)) = 0;

% Convert to matrices
sprime = cell2mat(MDP.info.sprime); % which state can be transitioned to
idx = sum(isnan(sprime),2) < size(sprime,2) == 1; % let non-nan states transition to themselves
sprime(idx,end) = find(idx);

T      = cell2mat(MDP.info.T);      % the probability of each sprime
R      = cell2mat(MDP.info.R);      % the reward at each sprime

% Run iterations
if makeplots
    figure
end
for i = 1:max_iter
    
    % Copy Q
    q = Q;

    % Precompute successor value
    switch policy
        case 'pessimism'
            V = w * nanmax(Q,[],2) + (1-w) * nanmin(Q,[],2);
        case 'softmax'
            V = nansum(Q .* (exp(beta*Q) ./ nansum(exp(beta*Q),2)),2);
        case 'min'
            V = nanmin(Q,[],2);
        case 'max'
            V = nanmax(Q,[],2);
        case 'mean'
            V = nanmean(Q,2);
    end
    
    % Correct the value of safe places
    if MDP.safety.on
        V(MDP.safety.loc_safe) = MDP.safety.sval;
    end
    
    % Assign V to each sprime
    V(end+1,:) = Inf;
    Vprime = sprime;
    Vprime(isnan(sprime)) = size(V,1);
    
    Vprime = V(Vprime);
    Vprime(Vprime==Inf) = NaN;
    
    % Update q value
    Q = R + gamma * T .* Vprime;
    
    % Compute delta
    delta = abs(Q - q);
    
    % Check for termination
    if all(delta < tol)
        break
    end
    
    MDP.V = V(1:end-1);
    
    % Draw progress
    if makeplots
        clf
        draw_map(MDP)
        drawnow
    end
    
end

% close(gcf)

if i == max_iter
    warning(['Q value computation reached maximum iterations (' num2str(max_iter) ')'])
end

MDP.Q = Q;

%% Solve for V

% get best sprime per state (according to max Q value)
max_idx = ones(size(Q,1),1)*5; % default to last action (transition to self)
for i = 1:size(Q,1)
    if any(~isnan(Q(i,:)))
        this_max = max(Q(i,:));
        if sum(Q(i,:)==this_max) == 1
            max_idx(i,1) = find(Q(i,:) == this_max);
        end
    end
end
% [~,max_idx] = nanmax(Q,[],2);

qpolicy = nan(MDP.nStates,1);
for st = 1:MDP.nStates
    qpolicy(st,1) = sprime(st,max_idx(st));
end

%% Compute policy

plan = MDP.start;

while true
   
    s = plan(end);
    
    if any(s == MDP.terminal)
        break
    end
    
    if any(plan == qpolicy(s))
        break;
    end
    
    plan = [plan qpolicy(s)];
    
end

MDP.plan = plan;

%% Plot

% Plot
if makeplots
    figure
    draw_map(MDP);
end

end