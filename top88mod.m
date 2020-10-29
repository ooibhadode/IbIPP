%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function xPhys = top88mod(nelx,nely,volfrac,penal,rmin,ft,F,fixeddofs,NonD,MusD,IM,q,beta,E0,v)
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
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
if ft == 1 || ft == 2
    xPhys = x; change_thresh = 1e-2;
else
    xTilde = x; change_thresh = 1e-4;
    xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
end
loop = 0; loopbeta = 0;
change = 1;
%% START ITERATION
while change > change_thresh && loop <= 300
  loop = loop + 1; loopbeta = loopbeta+1;
  %% FE-ANALYSIS
  if strcmpi(IM,'SIMP')
      sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  else
      sK = reshape(KE(:)*(xPhys(:)'./(1+q*(1-xPhys(:)')))*E0,64*nelx*nely,1);
  end
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS        % changed
  c = 0; dc = 0;
  for i = 1:size(F,2)
      Ui = U(:,i);
      ce=reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
      if strcmpi(IM,'SIMP')
          c = c + sum((Emin+xPhys(:)'.^penal*(E0-Emin))*ce(:));
          dc = dc-penal.*xPhys.^(penal-1).*(E0-Emin).*ce;
      else
          c = c + sum((xPhys(:)'./(1+q*(1-xPhys(:)')))*(E0)*ce(:));
          dc = dc-(((1+q*(1-xPhys))+q*xPhys)./(1+q*(1-xPhys)).^2).*E0.*ce;
      end
  end
  dc = reshape(dc,nely,nelx);
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  elseif ft == 3
    dxHv = beta*exp(-beta*xTilde)+exp(-beta);
    dc(:) = H*((dc(:).*dxHv(:))./Hs);
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    elseif ft == 3
      xTilde = H*(xnew(:)./Hs);
      xPhys(:) = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    end
    xPhys = xPhys(:); xPhys(MusD) = 1; xPhys(NonD) = 0; xPhys = reshape(xPhys,nely,nelx);  % New line
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    sum(xPhys(:))/(nelx*nely-length(NonD)),change);
  %% PLOT DENSITIES
  figure (1)
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
 if ft == 3 && beta < 512 && (mod(loopbeta,50)==0 || change <=1E-2)
     beta = 2*beta;
 else
 end
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
