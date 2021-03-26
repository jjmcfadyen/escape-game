function [outcome] = mdp_dynamic(MDP,dyn)

MDP.frameRate = dyn.frameRate;

% Create controllers
agent = controller(MDP,'agent','initiate');
predator = controller(MDP,'predator','initiate');

% Update dynamic variables
agent.speed     = dyn.agentSpeed;
predator.speed  = dyn.predatorSpeed;
agent.lookahead = dyn.lookahead;

% Determine max number of frames for this trial
frames = dyn.timeLimit*dyn.frameRate;

%% Run through frames

figure
set(gcf,'position',[2258 -331 868 647])

f = 0;
while f <= frames
   
    f = f + 1;
    
    % swap information between agent & predator
    predator.agent_hiding = agent.hiding;
    predator.agent(1)     = agent.start;
    agent.predator(1)     = predator.start;
    
    % check if died/won
    if predator.start == agent.start || agent.start == agent.reward(1)
        if predator.start == agent.start
            outcome = 'died';
        elseif agent.start == agent.reward(1)
            outcome = 'won';
        end
        break
    end
    
    % evaluate agent's environment
%     if isempty(agent.moving)
        agent = controller(agent,'agent','evaluate');
%     end
%     if isempty(predator.moving)
        if agent.hiding
            predator.plan = predator.plan(find(predator.plan==predator.start):end,:);
            if predator.plan(end) == predator.start
                predator = controller(predator,'predator','evaluate');
            end
        else
            predator = controller(predator,'predator','evaluate');
        end
%     end
    
    % set movement
    if isempty(agent.moving) % if the agent has no trajectory...
        agent = controller(agent,'agent','setmovement');
    end
    if isempty(predator.moving) % if the agent has no trajectory...
        predator = controller(predator,'predator','setmovement');
    end
    
    % move agent
    if ~isempty(agent.moving)
        agent.position = agent.moving(1,:);
        agent.moving = agent.moving(2:end,:);
        if isempty(agent.moving)
            agent.start = agent.map(agent.position(1),agent.position(2));
            agent.hiding = any(agent.start == agent.safety.loc_safe);
        end
    end
    
    % move predator
    if ~isempty(predator.moving)
        predator.position = predator.moving(1,:);
        predator.moving = predator.moving(2:end,:);
        if isempty(predator.moving)
            predator.start = predator.map(predator.position(1),predator.position(2));
        end
    end
    
    % plot
    clf;
    draw_map(agent);
    scatter(agent.position(1,2),agent.position(1,1),...
        100,'markerfacecolor','b','markeredgecolor','k'); hold on
    scatter(predator.position(1,2),predator.position(1,1),...
        100,'markerfacecolor','r','markeredgecolor','k'); hold on
    caxis([-11 10])
    if agent.hiding
        title('hiding')
    else
        title('moving')
    end
    drawnow

end


end