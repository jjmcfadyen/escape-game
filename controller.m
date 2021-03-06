function Y = controller(X,type,job)

MDP = X;

switch job
    
    case 'initiate'

        %% Initiate variables for the agent OR predator
        
        % Copy X
        Y = X;
        
        switch type
            case 'agent'
                
                % Set variables
                Y.cansee    = [];
                Y.plan      = [];
                Y.predator  = X.punishment;
                Y.hiding    = false;
                Y.moving    = [];
                Y.cone      = Y.agent_cone;
                
                % Remove fields
                if isfield(Y,'V')
                    Y = rmfield(Y,'plan');
                end
                if isfield(Y,'plan')
                    Y = rmfield(Y,'plan');
                end
                
            case 'predator'
                
                % Set variables
                Y.start     = Y.punishment(1);  
                Y.agent     = [X.start X.reward(2)];
                Y.plan      = [];
                Y.moving    = [];
                Y.cone      = Y.predator_cone;
                Y.terminal  = [];
                
                % Remove fields
                remove_fields = {'punishment','reward','V','plan'};
                for f = 1:length(remove_fields)
                    if isfield(Y,remove_fields{f})
                        Y = rmfield(Y,remove_fields{f});
                    end
                end
        end
        
        % Set common features
        Y.position  = u_coordswitch(Y.start,Y);
        
    case 'evaluate'
        
        %% Evaluate the map
        
        switch type
            case 'agent'
                
                % Check what can currently be seen
                if X.cone ~= Inf
                    [~,labels]  = u_viscone(X.start,X);
                    X.cansee    = labels;
                else
                    X.cansee    = 1:X.nStates;
                end
                
                % Can the predator be seen?
                if any(X.cansee == X.predator(1)) % yes
                    X.punishment = X.predator;
                    X.terminal      = [X.reward(1) X.punishment(1)];
                elseif isfield(X,'punishment') % no
                    X.terminal      = setdiff(X.terminal,X.punishment(:,1));
                    X               = rmfield(X,'punishment');
                end
                
                % If multiple rewards, focus on the first one
                if size(X.reward,1) > 1
                    X.reward = X.reward(1,:);
                end
                
                % Update entire value map (but only if different from before)
                X = mdp_build(X);
                
                if ~isfield(X,'V') || u_valdiff(X,MDP)
                    X = mdp_valueIteration(X,false);
                else
                    pos = u_coordswitch(round(X.position),X); 
                    if pos ~= X.plan(1) % update position in plan
                        X.plan = X.plan(find(X.plan==pos):end);
                    end
                end
                
                % Check if hiding
                if X.safety.on
                    if length(X.plan) > 2
                        nextmove = X.plan(2);
                    else
                        try
                            nextmove = X.plan(1);
                        catch
                            X = mdp_valueIteration(X,false);
                            nextmove = X.plan(1);
                        end
                    end
                    if X.start == nextmove && any(X.start == X.safety.loc_safe)
                        X.hiding = true;
                    end
                end

            case 'predator'

                % Check what can currently be seen
                if X.cone ~= Inf
                    [~,labels]  = u_viscone(X.start,X);
                    X.cansee    = labels;
                else
                    X.cansee    = 1:X.nStates;
                end
                
                % Can the agent be seen?
                if any(X.cansee == X.agent(1)) &&  ~X.agent_hiding % yes
                    
                    % Compute path to agent
                    X.plan = mdp_astar(X.start,X.agent(1),X);
                    
                else
                    
                    tmp = X;
                    if X.agent_hiding % Make their safe location unavailable
                        tmp.map(tmp.map==X.agent(1)) = NaN;
                        tmp = mdp_build(tmp);
                    end
                    
                    % Compute path to random location   
                    loc_rand = tmp.viableStates(randi(tmp.nViable));
                    X.plan = mdp_astar(tmp.start,loc_rand,tmp);
                    
                end
                
                X.plan = [X.start; X.plan];
                
        end
        
        % Output
        Y = X;
        
    case 'setmovement'
    
        % get current position
        position = X.position;
        
        % extract next planned move
        if length(X.plan) > 1
            destination = u_coordswitch(X.plan(2),X);
        else
            destination = position;
        end
        
        % see how many frames until they will reach next square
        move_frames = 0:0+(X.frameRate/X.speed)-1;
        trajectory = [linspace(position(1),destination(1),length(move_frames))',...
                      linspace(position(2),destination(2),length(move_frames))'];
        
        % update trajectory
        X.moving = trajectory;
        
        % output
        Y = X;
        
    case 'checkstatus'
        
        switch type
            case 'agent'
                % Check if safe
                if any(ismember(Y.safety.loc_safe,Y.start))
                    Y.issafe = true;
                else
                    Y.issafe = false;
                end
        end
        
end

end