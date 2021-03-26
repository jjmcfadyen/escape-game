function [coords,labels] = u_viscone(A,MDP)

    % check coordinates
    if length(A) == 1 % conver state label to [x,y] coordinate
         A = u_coordswitch(A,MDP);
    end

    % get viewing window
    rows = A(1)-MDP.cone:A(1)+MDP.cone;
    cols = A(2)-MDP.cone:A(2)+MDP.cone;
    
    rows(rows > MDP.dim(1) | rows <= 0) = [];
    cols(cols > MDP.dim(2) | cols <= 0) = [];
    
    window = MDP.map(rows,cols);
    [I,J] = find(window == MDP.map(A(1),A(2)));
    A = [I,J];
    
    window(~isnan(window)) = 0;
    window(isnan(window)) = 1;
    
    % find maximum distance
    maxdist = max([norm(A), norm(A-[1 size(window,2)]), norm(A-[size(window,1) 1]), norm(A-numel(window))]);
    angles  = 1:360; % use smaller increment to increase resolution
    
    % generate set of points uniformly distributed around start point
    endpoints = bsxfun(@plus, maxdist*[cosd(angles)' sind(angles)'], A);

    intersec = [];
    for k = 1:numel(angles)
        [CX,CY,C] = improfile(window,[A(2), endpoints(k,1)],[A(1), endpoints(k,2)]);
        idx = find(C);
        intersec(k,:) = [CX(idx(1)), CY(idx(1))];
    end
    
    fov = roipoly(window,intersec(:,1),intersec(:,2));
    fov(window==1) = 0;
    
    % retrieve state positions
    labels = MDP.map(rows,cols);
    labels = labels(fov);
    
    coords = nan(length(labels),2);
    for h = 1:length(labels)
        [I,J] = find(MDP.map == labels(h));
        coords(h,:) = [I,J];
    end
    
end