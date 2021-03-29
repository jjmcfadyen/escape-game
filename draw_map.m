function draw_map(MDP,opts)


%% Plot

% Set up figure
set(gca,'ydir','reverse')
xlim([0.5 MDP.dim(1)+0.5])
ylim([0.5 MDP.dim(2)+0.5])

title('Grid World')
set(gca,'xtick',1:MDP.dim(1))
set(gca,'ytick',1:MDP.dim(2))

% Draw value map (if applicable)
if isfield(MDP,'V')
    
    % Reshape values to match map
    V = reshape(MDP.V,MDP.dim(1),MDP.dim(2))';
    V(isnan(MDP.map)) = NaN;

    % Fill in terminal states
%     if MDP.safety.on && ~isempty(MDP.safety.found)
%         for i = 1:length(MDP.safety.found)
%             [I,J] = find(MDP.map == MDP.safety.found(i));
%             V(I,J) =MDP.safety.sval;
%         end
%     end
    
    if isfield(MDP,'reward')
        for i = 1:size(MDP.reward,1)
            [I,J] = find(MDP.map == MDP.reward(i,1));
            V(I,J) = MDP.reward(i,2);
        end
    end
    
    if isfield(MDP,'punishment')
        [I,J] = find(MDP.map == MDP.punishment(1));
        V(I,J) = MDP.punishment(2);
    end

    % set colormap
    imagesc(V); hold on
    
    cmap = [0 0 0; u_colours(1000,'viridis')];
    colormap(cmap);

end

% Draw gridlines
% for r = 0.5:1:MDP.dim(1)+0.5
%     plot([r r],[0.5 MDP.dim(2)+0.5],'k','linewidth',.5); hold on
%     plot([0.5 MDP.dim(2)+0.5],[r r],'k','linewidth',.5); hold on
% end

% Put border around terminal states
for t = 1:length(MDP.terminal)
    [row,col] = find(MDP.map == MDP.terminal(t));
    plot(repmat(col-.5,2,1),[row-.5 row+.5],'w','linewidth',3);
    plot(repmat(col+.5,2,1),[row-.5 row+.5],'w','linewidth',3);
    plot([col-.5 col+.5],repmat(row-.5,2,1),'w','linewidth',3);
    plot([col-.5 col+.5],repmat(row+.5,2,1),'w','linewidth',3);
end

% Show policy
if isfield(MDP,'plan')
    for t = 1:length(MDP.plan)-1
        p1 = [0 0];
        p2 = [0 0];
        [p1(2) p1(1)] = find(MDP.map == MDP.plan(t));
        [p2(2) p2(1)] = find(MDP.map == MDP.plan(t+1));
        dp = (p2-p1) * .75;
        tmp = quiver(p1(1),p1(2),dp(1),dp(2),0,'w','maxheadsize',50);
    end
end

% Show safe locations
if MDP.safety.on
    if isfield(MDP.safety,'found')
        for i = 1:length(MDP.safety.found)
            [I,J] = find(MDP.map == MDP.safety.found(i));
            scatter(J,I,100,'markerfacecolor','y','markeredgecolor','k');
        end 
    end
end

% Tidy
set(gca,'ticklength',[0 0])

end