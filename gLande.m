function [gLande] = gLande(L,S)

if length(L) ~= length(S)
    disp('gLande(): L & S must be equal length')
    return
end
J = L + S;
gLande = (3/2)+((S.*(S+1)-L.*(L+1))./(2.*J.*(J+1)));
end