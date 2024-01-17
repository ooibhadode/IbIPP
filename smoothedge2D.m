function [vxPhys,xg,top,tol] = smoothedge2D(vxPhys,Hn,Hns,nelx,nely,nele,nnele,nodex,nodey,fnx,fny,beta,ngrid)

xn=reshape((Hn*vxPhys(:)./Hns),nely+1,nelx+1);
%% UPDATE POINT DESNIGY BY A HEAVISIDE SMOOTH/STEP FUNCTION
xg=interp2(nodex,nodey,xn,fnx,fny,'linear'); clear xn
l1=0; l2=1;
while (l2-l1)>1e-5
ls=(l1+l2)/2;
xgnew=max(0.001,(tanh(beta*ls)+tanh(beta*(xg-ls)))/(tanh(beta*ls)+tanh(beta*(1-ls))));
if sum(sum(xgnew))/((ngrid*nelx+1)*(ngrid*nely+1))-sum(vxPhys(:))/(nele)>0, l1=ls; else, l2=ls; end
end
clear l1 l2
%% ASSEMBLE GRID POINTS TO ELEMENTS
vxPhys(:)=0;
Terr=0; Tm=[];
for i=1:nelx
    for j=1:nely
        e=(i-1)*nely+j;
            for i1=ngrid*(i-1)+1:ngrid*i+1
                for j1 = ngrid*(j-1)+1:ngrid*j+1
                    Tm = [Tm;xgnew(j1,i1)];
                    vxPhys(e)=vxPhys(e)+xgnew(j1,i1);
                end
            end
        if min(Tm)>0.001 && max(Tm)<1, Terr=Terr+1; end
        Tm=[];        
    end
end
vxPhys=vxPhys/(ngrid+1)^2;
tol=Terr/nnele;
top = xg-ls;