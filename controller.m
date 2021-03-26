function Y = controller(X,type,job)

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
                
                X.plan = [];
                
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
                    if X.lookahead > 0
                        predator_path = mdp_astar(X.predator(1),round(X.position),X);
                        if X.lookahead < length(predator_path)
                            X.punishment(1) = predator_path(X.lookahead);
                        end
                    end
                    X.terminal      = [X.reward(1) X.punishment(1)];
                elseif isfield(X,'punishment') % no
                    X               = rmfield(X,'punishment');
                    X.terminal      = X.reward(1);
                end
                
                % Can any safe locations be seen?
                if X.safety.on
                    found = labels(ismember(labels,X.safety.loc_safe));
                    if ~isempty(found)
                        X.safety.found = unique([X.safety.found found]);
                    end
                end
                
                % Update entire value map
                X = mdp_build(X);
                X = mdp_valueIteration(X,false);
                
                % Check if hiding
                if X.safety.on
                    if length(X.plan) > 2
                        nextmove = X.plan(2);
                    else
                        nextmove = X.plan(1);
                    end
                    if X.start == nextmove && any(X.start == X.safety.found)
                        X.hiding = true;
                    end
                end

            case 'predator'

                % Check if agent is hiding or not
                if ~X.agent_hiding % no
                    
                    % Compute path to agent
                    X.plan = mdp_astar(X.start,X.agent(1),X);
                    
                else % yes
                    
                    % Make their safe location unavailable
                    tmp = X;
                    tmp.map(tmp.map==X.agent(1)) = NaN;
                    tmp = mdp_build(tmp);
                    
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