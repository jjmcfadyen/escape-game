clear all
close all
clc

%% Define parameters

gridworld   = zeros(5,5);

% Successor representations
Mepsilon     = 0;             % the transition uncertainty between states (max is 0.5; anything higher, they will do the OPPOSITE of what is likely)
policy      = 'pessimism';   % policy        (pessimism, softmax, min, max, mean)
gamma       = 0.98;          % gamma         (temporal discounting factor)
beta        = 5;             % beta          (inverse temperature - softmax only)
w           = 0.7;           % reward bias   (pessimism only)
tol         = 1e-4;
max_iter    = 100;

% Visibility setting
cone        = 3;    % how many squares ahead the agent can see

% Reward settings
reward              = [];
reward.val          = 10;   
reward.num          = 2;         % how many switches ( + exit)

% Safety settings
safety              = [];
safety.on           = true;      % whether safe locations are available or not
safety.sdensity     = .05;       % proportion of the environment with safe locations (0 to 1)
safety.sval         = 0;         % perceived value of safety zone
safety.found        = [];        % log of safety locations found
safety.method       = 'merged';  % how safety locations are dealt with ('merged')

% Dynamic settings 
dyn                  = [];
dyn.on               = true;    % enables you to produce a score (escape or death)
dyn.timeLimit        = 60;      % how long the agent has to escape
dyn.predatorSpeed    = 1.5;      % how fast the predator can move (in moves per second)
dyn.agentSpeed       = 1;       % how fast the agent can move (in moves per second)
dyn.predatorFunction = 'astar'; % 'astar' or 'softmax'
dyn.frameRate        = 60;      % update rate of movement
dyn.lookahead        = true;    % does the agent imagine deaths along future policies, given predator speed?

%% Initiate MDP

MDP = [];

MDP.nStates     = numel(gridworld);
MDP.dim         = size(gridworld);
MDP.map         = reshape(1:MDP.nStates,MDP.dim(1),MDP.dim(2))';
MDP.reward      = reward;
MDP.safety      = safety;
MDP.cone        = cone;

MDP = u_maze(MDP,dyn);

MDP.reward          = [MDP.reward, repmat(10,length(MDP.reward),1)];       % set value of terminal reward state
MDP.punishment(2)   = -10;              % set value of terminal predator state
MDP.terminal        = [MDP.reward(end,1) MDP.punishment(1)];  % add terminal state as predator state

% Log other settings
MDP.epsilon     = epsilon;
MDP.policy      = policy;
MDP.gamma       = gamma;
MDP.beta        = beta;
MDP.w           = w;
MDP.tol         = tol;
MDP.max_iter    = max_iter;

MDP = mdp_build(MDP);

%% Dangerous environment

if dyn.on
    
%     % randomly assign predator location
%     MDP = u_maze(MDP,'predator'); % 'all' or 'predator'
    
    % run the chase scene
    MDP.w = 0.8;
    [nonanx_outcome] = mdp_dynamic(MDP,dyn);
    
    MDP.w = 0.2;
    [anx_outcome] = mdp_dynamic(MDP,dyn);
    
else
    MDP = mdp_valueIteration(MDP,true); % build a static map from the beginning state
end
