function output = u_valdiff(A,B)
% Tells you if the transition or reward matrices are different between 2
% MDPs

AT = cell2mat(A.info.T);
AR = cell2mat(A.info.R);

BT = cell2mat(B.info.T);
BR = cell2mat(B.info.R);

AT(isnan(AT)) = Inf;
AR(isnan(AR)) = Inf;
BT(isnan(BT)) = Inf;
BR(isnan(BR)) = Inf;

output = any(AT(:) ~= BT(:)) || any(AR(:) ~= BR(:));

end