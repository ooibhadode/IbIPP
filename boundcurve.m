function [BN1,ang] = boundcurve(edofMat1,nely,nelx,xp)
%% Obtain boundary nodes
b_n = zeros(1,(nelx+1)*(nely+1));
for i = 1:(nelx+1)*(nely+1)
    [row, ~] = find(edofMat1 == i); mean_x = mean(xp(row));
    if mean_x < 1 && mean_x > 0
        b_n(i) = i;
    else
    end
end
bn1 = b_n(:); bn1(b_n == 0)=[]; BN1 = bn1;
sx = reshape(repmat(1:nelx+1,nely+1,1),(nely+1)*(nelx+1),1);
sy = reshape(repmat((nely+1:-1:1)',1,nelx+1),(nely+1)*(nelx+1),1);
dist = zeros(length(BN1),length(BN1));
bn_x = sx(BN1)'; bn_y = sy(BN1)'; positn = zeros(length(BN1),2);

%% Obtain corresponding nodal angles
ang = zeros(length(BN1),1);
for i = 1:length(BN1)
    dist(i,:) = sqrt((bn_x(i)-bn_x).^2+(bn_y(i)-bn_y).^2);
    distt = dist(i,:); dis3t = mink(distt,3);
    posit = find(distt == dis3t(2));
    if length(posit) == 1
      positn(i,1:2) = [posit find(distt == dis3t(3),1,'first')];
    elseif length(posit) > 2
      positn(i,1:2) = posit(1:2);
    else
      positn(i,1:2) = posit;
    end
    ang(i) = atand((bn_y(positn(i,1))-bn_y(positn(i,2)))/(bn_x(positn(i,1))...
        -bn_x(positn(i,2))));
end
