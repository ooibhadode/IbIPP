%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function [xPhys,xg] = semdot(nelx,nely,volfrac,rmin,F,fixeddofs,NonD,MusD,beta,E0,v,ER,ngrid)
%% MATERIAL PROPERTIES
E0 = E0*1;
Emin = E0*(1e-9);
nu = v;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);                  % Deleted
U = zeros(2*(nely+1)*(nelx+1),size(F,2));                   % Changed
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); % Deleted
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
[Hn,Hns]=HnHns(nelx,nely,1);
%% ELEMENTAL NODES AND COORDINATES
[nodex,nodey]=meshgrid(0:nelx,0:nely);
[fnx,fny]=meshgrid(0:1/ngrid:nelx,0:1/ngrid:nely);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x; change_thresh = 1e-2;
loop = 0; loopbeta = 0;
change = 1; tol = 1; tolx = 1e-3; tol_thresh = 1e-3; nloop = 300;
%% START ITERATION
while (change > tolx || tol>tol_thresh) && loop < nloop
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(xPhys(:)'*E0+(1-xPhys(:))'*Emin*E0),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS        % changed
  c = 0; dc = 0;
  for i = 1:size(F,2)
      Ui = U(:,i);
      ce=reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
      c = c + sum(((1-xPhys(:))*Emin+xPhys(:)).*E0.*ce(:));
      dc = dc-((1-xPhys)*Emin+xPhys).*E0.*ce;
  end
  dc = reshape(dc,nely,nelx);
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  dc(:) = H*(dc(:)./Hs);
  dv(:) = H*(dv(:)./Hs);
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    xPhys(:) = (H*xnew(:))./Hs;
    xPhys = xPhys(:); xPhys(MusD) = 1; xPhys(NonD) = 0; xPhys = reshape(xPhys,nely,nelx);  % New line
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew; nele = nelx*nely; nnele = nelx*nely-length(NonD);
  [xPhys,xg,top,tol] = smoothedge2D(xPhys,Hn,Hns,nelx,nely,nele,nnele,nodex,nodey,fnx,fny,beta,ngrid);
  %% PRINT RESULTS
   disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(xPhys(:))/nnele) ...
        ' Ch.: ' sprintf('%6.3f',change) ...
        ' Tol.: ' sprintf('%6.3f',tol)]);
  %% PLOT DENSITIES
  figure (1)
  contourf(fnx, flipud(fny), top, [0 0]); axis equal tight off; drawnow
  beta = beta+ER;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
