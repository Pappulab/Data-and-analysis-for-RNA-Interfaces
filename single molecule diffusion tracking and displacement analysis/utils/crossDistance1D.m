function [Vind, d] = crossDistance1D(U,V,threshold)
% U, V must be i-by-2 and j-by-2 x-y coordinated locations
if ~isempty(U) && ~isempty(V)
    Ux = U(:,1); Uy = U(:,2);
    Vx = V(:,1); Vy = V(:,2);
    [UX,VX] = meshgrid(Ux,Vx);
    [UY,VY] = meshgrid(Uy,Vy);
    dMat = sqrt((UX-VX).^2+(UY-VY).^2);
    [d,Vind] = min(dMat,[],'all');
    if d>threshold
        Vind = [];
        d = [];
    end
else
    Vind = [];
    d = [];
end
end