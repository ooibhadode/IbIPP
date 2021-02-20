function datatostl(nelx,nely,xPhy,tx,form,ht,lr,theta,symm)
xp = xPhy; xp(xp >= 0.5) = 1; 
%% Check for symmetry requirements
if strcmpi(symm,'left')
    xp(:,sum(xp,1)<0.01*nely) = [];
    xp = [fliplr(xp) xp]; 
    nelx = size(xp,2); nely = size(xp,1);
elseif strcmpi(symm,'bottom')
    xp(sum(xp,2)<0.01*nelx,:) = [];
    xp = [xp; flipud(xp)];
    nelx = size(xp,2); nely = size(xp,1);
elseif strcmpi(symm,'right')
    xp(:,sum(xp,1)<0.01*nely) = [];
    xp = [xp fliplr(xp)];
    nelx = size(xp,2); nely = size(xp,1);
elseif strcmpi(symm,'top')
    xp(sum(xp,2)<0.01*nelx,:) = [];
    xp = [flipud(xp); xp];
    nelx = size(xp,2); nely = size(xp,1);
else 
end 
%% Extrude or revolve
xp(:,sum(xp,1)<0.01*nely) = []; xp(sum(xp,2)<0.01*nelx,:) = [];
[nely,nelx] = size(xp);
xp = [zeros(1,nelx);xp;zeros(1,nelx)]; xp = [zeros(nely+2,1),xp,zeros(nely+2,1)];
if strcmpi(form,'extrude')
    nelz = round(ht*min(nely,nelx));
    df = zeros(nely+2,nelx+2,nelz+2); df(:,:,2:nelz+1) = repmat(xp,[1 1 nelz]);
elseif strcmpi(form,'revolve')
    df1 = [zeros(size(xp,1),round(lr*nelx)) xp];
    if theta == 360
        df1 = xp(:,2:end);
    end 
    df = revolve2D(df1,theta);
end    
fv = isosurface(df,0.5); slt = split(tx,'.'); mkdir(slt{1}); stlwrite([slt{1} '\' tx],fv);