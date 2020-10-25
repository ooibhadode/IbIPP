function [F,fixeddofs,NonD,MusD,volfrac,edofMat1] = loadandsupport(nelx,...
    nely,Fmag,Fang,Pmag,dom,pre_support,pre_load,vol_frac)
nele = nelx*nely;                                                          % Number of elements in domain
%% Prepare mesh
nodenr = reshape(1:(1+nelx)*(1+nely),(1+nely),(1+nelx));                   % node numbers in domain
nodenrs1 = nodenr-1; nodenrs1(end,:) = []; nodenrs1(:,end) = [];
edofMat1 = nodenrs1(:) + [2 2+(nely+1) 2+nely 1];                          % connectivity matrix for node numbers
nodenrs = 2*(nodenr-1); nodenrs(end,:) = []; nodenrs(:,end) = [];
edofMat = nodenrs(:) + [3 4 5+2*nely 6+2*nely 3+2*nely 4+2*nely 1 2];      % connectivity matrix for degrees of freedom

%% Load definition 
dom1 = dom(:); 
% Forces
Fp = zeros(2*(nelx+1)*(nely+1),length(Fmag));                              % force matrix initialization
for i = 1:length(Fmag)
    F_ele = find(dom1 == 20+i-1);                                          % obtain the force elements
    if ~isempty(F_ele)
        fnodes11 = edofMat1(F_ele,:);                                      % obtain elements' nodes 
        fnodes11 = unique(fnodes11(:));
        fnodes1 = fnodes11(round(length(fnodes11)/2));                     % obtain node at the ~centre of elements
        [rw,~] = find(edofMat1 == fnodes1);
        if mean(dom1(rw)) < dom1(F_ele(1))                                 % check that the node is at the edge
            fn1 = fnodes11(fnodes11 < fnodes11(1)+nely+1);                  
            fnodes1 = fnodes1-floor(0.5*length(fn1));                      % Ensure it is centered
        end 
        Fp(2*fnodes1-1,i) = Fmag(i)*sind(Fang(i));                         % resolve x-force component
        Fp(2*fnodes1,i) = Fmag(i)*cosd(Fang(i));                           % resolve y-force component
    else 
    end 
end

% Pressure
P = zeros(2*(nelx+1)*(nely+1),length(Pmag));                               % Initialize pressure matrix
dom2 = dom; dom2(dom2 > 1) = 1; dom2 = dom2(:);                            % prepare domain for boundary ID

if ~isempty(Pmag)
    [BN1,ang] = boundcurve(edofMat1,nely,nelx,dom2);                       % Boundary ID to obtain boundary nodes and corresponding angles
    for i = 1:length(Pmag)
        P_ele = find(dom1 == 30+i-1);                                      % obtain pressure elements
        if ~isempty(P_ele)
            pnodes1 = edofMat1(P_ele,:); pnodes1 = unique(pnodes1(:));
            pnodes = intersect(pnodes1,BN1);                               % obtain pressure nodes
            Pmag1 = Pmag(i)/length(pnodes);                                % obtain nodal equivalent force
            for j = 1:length(pnodes)
                [row1, ~] = find(edofMat1 == pnodes(j)+nely);              % obtain the mean density of elements surrounding a node 
                mean_x1 = mean(dom2(row1));                                %that is nely in front of pnodes(j)
                [row2, ~] = find(edofMat1 == pnodes(j)-(nely+2));          % similar to previous comment
                mean_x2 = mean(dom2(row2));
                [row3, ~] = find(edofMat1 == pnodes(j)-(nely+1));          % similar to previous comment
                mean_x3 = mean(dom2(row3));
                [row4, ~] = find(edofMat1 == pnodes(j)-1);                 % similar to previous comment
                mean_x4 = mean(dom2(row4));
                Ang = ang(BN1 == pnodes(j));                               % obtain angle of pressure node under investigation
                if Ang<0 && mean_x1>0
                    P(2*pnodes(j)-1,i) = Pmag1*cosd(abs(Ang));             % resolve equivalent force to its components
                    P(2*pnodes(j),i) = Pmag1*sind(abs(Ang));
                elseif Ang<0 && mean_x1==0
                    P(2*pnodes(j)-1,i) = -Pmag1*cosd(abs(Ang));            % resolve equivalent force to its components
                    P(2*pnodes(j),i) = -Pmag1*sind(abs(Ang));
                elseif Ang>0 && mean_x2>0
                    P(2*pnodes(j)-1,i) = -Pmag1*cosd(abs(Ang));            % resolve equivalent force to its components
                    P(2*pnodes(j),i) = Pmag1*sind(abs(Ang));
                elseif Ang>0 && mean_x2 == 0
                    P(2*pnodes(j)-1,i) = Pmag1*cosd(abs(Ang));             % resolve equivalent force to its components
                    P(2*pnodes(j),i) = -Pmag1*sind(abs(Ang));
                elseif abs(Ang) == 90 && mean_x3>0
                    P(2*pnodes(j)-1,i) = -Pmag1;                           % resolve equivalent force to its components
                elseif abs(Ang) == 90 && mean_x3 == 0
                    P(2*pnodes(j)-1,i) = Pmag1;                            % resolve equivalent force to its components
                elseif abs(Ang) == 0 && mean_x4>0
                    P(2*pnodes(j),i) = Pmag1;                              % resolve equivalent force to its components
                elseif abs(Ang) == 0 && mean_x4 == 0
                    P(2*pnodes(j),i) = -Pmag1;                             % resolve equivalent force to its components
                else
                end 
            end 
        end 
    end 
end 
F = [Fp,P]; F(:,all(F == 0)) = [];                                          % concatenate foroce and pressure matrices, delete all-zero columns         

%% Support definition
% Fixeddofs and freedofs
fixeddofs1 = find(dom1 == 51); fixeddofs2 = find(dom1 == 52);               % obtain elements with fixed DOFs
fixeddofs3 = find(dom1 == 53); 
if ~isempty(fixeddofs1)
    fixeddofs_1 = edofMat(fixeddofs1,:); fixeddofs_1 = unique(fixeddofs_1(:)); % elements with all fixed DOFs
else 
    fixeddofs_1 = [];
end 
if ~isempty(fixeddofs2)
    fixeddofs_2 = edofMat1(fixeddofs2,:); fixeddofs_2 = unique(fixeddofs_2(:)); % elements with fixed x-DOFs
    fixeddofs_2 = 2*fixeddofs_2-1;
else
    fixeddofs_2 = [];
end 
if ~isempty(fixeddofs3)
    fixeddofs_3 = edofMat1(fixeddofs3,:); fixeddofs_3 = unique(fixeddofs_3(:)); % elements with fixed y-DOFs
    fixeddofs_3 = 2*fixeddofs_3;
else
    fixeddofs_3 = [];
end 
fixeddofs = union(fixeddofs_1,union(fixeddofs_2,fixeddofs_3));             % compute a collective fixeddofs

%% Preserved elements `
if pre_support == 0
    MusD = find(dom1 == 4);                                                % preserved elements
elseif pre_support == 1
    MusD = union(find(dom1 == 4),find(dom1 == 51));                        % preserved and fixed support elements
elseif pre_support == 2
    MusD = union(find(dom1 == 4),find(dom1 == 52));                        % preserved and x-fixed support elements
elseif pre_support == 3
    MusD = union(find(dom1 == 4),find(dom1 == 53));                        % preserved and y-fixed support elements
elseif pre_support == 4
    MusD = union(find(dom1 == 4),union(find(dom1 == 51),find(dom1 == 52)));% preserved, fixed and x-fixed support elements
elseif pre_support == 5
    MusD = union(find(dom1 == 4),union(find(dom1 == 51),find(dom1 == 53)));% preserved, fixed and y-fixed support elements
elseif pre_support == 6
    MusD = union(find(dom1 == 4),union(find(dom1 == 52),find(dom1 == 53)));% preserved, x-fixed and y-fixed support elements
elseif pre_support == 7 
    MusD = union(find(dom1 == 4),union(find(dom1 == 51),union(find(dom1 ...% preserved, fixed, x-fixed and y-fixed support elements
        == 52),find(dom1 == 53))));
else 
    MusD = [];
end
% To preserve load elements
if pre_load == 1
    MusD = union(MusD,intersect(find(dom1 >= 20),find(dom1 < 30)));        % declare preserve force elements
elseif pre_load == 2
    MusD = union(MusD,intersect(find(dom1 >= 30),find(dom1 < 40)));        % declare preserved pressure elements
elseif pre_load == 3
    MusD = union(MusD,union(intersect(find(dom1 >= 20),find(dom1 < 30)),...
        intersect(find(dom1 >= 30),find(dom1 < 40))));                     % declare preserved force & pressure elements
else 
end 
%% Non design elements
NonD = find(dom1 == 0);                                                    % declare non-design domain elements

%% New volume fraction
volfrac = vol_frac*(nele-length(NonD))/nele;                               % compute new volume fraction