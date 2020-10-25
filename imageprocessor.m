function [dom,nely] = imageprocessor(domain,nelx)
A = imread(domain);                                                         % Read image file
[rows,cols,~] = size(A); nely = round((rows/cols)*nelx);                    % obtain image resolution and nely
top = cell(nely,nelx); dom = zeros(nely,nelx);                              % initialize design domain with nely by nelx size  
for i = 1:nely
    for j = 1:nelx
        top{i,j} = A(round((i/(nely+1))*rows):round(((i+1)/(nely+1))*rows),... % obtain image area for every element in discretization
            round((j/(nelx+1))*cols):round(((j+1)/(nelx+1))*cols),:);
        top1 = top{i,j};
        rc = mean(mean(top1(:,:,1))); gc = mean(mean(top1(:,:,2)));         % obtain RGB value for each elemental area
        bc = mean(mean(top1(:,:,3)));
        if [rc,gc,bc] < 200                                                 % Black = element is in design domain
            dom(i,j) = 1;
        elseif rc >= 200 && gc < 5 && bc < 5                                % Red = force carrying element
            dom(i,j) = 20+floor((rc-200)/5);
        elseif rc >= 200 && rc <= 230 && gc >= 100 && gc <= 150 && bc < 5   % Orange = pressure carrying element 
            dom(i,j) = 30+floor((gc-100)/5);
        elseif rc <= 100 && gc >= 200 && bc <= 100                          % Green = preserved element
            dom(i,j) = 4;
        elseif rc < 100 && gc <= 100 && bc >= 200                           % Blue = element with fixed nodes
            dom(i,j) = 51;
        elseif rc >= 100 && rc <= 150 && gc < 100 && bc >= 200             % Purple = element with nodes fixed in x
            dom(i,j) = 52;
        elseif rc < 100 && gc >= 200 && bc >= 200                           % Cyan = element with nodes fixed in y
            dom(i,j) = 53;
        else                                                                % White = element is in non-design domain
        end 
    end
end