function Y = u_coordswitch(X,MDP)

% Convert state number to XY coordinate
if length(X) == 1 
    [I,J] = find(MDP.map == X);
    Y = [I,J];
% Convert XY coordinate to state number
elseif length(X) == 2 
    Y = MDP.map(X(1),X(2));
else
    error('This is neither a state number nor an XY coordinate')
end

end