%%%%%%%%%%%%%%%%%% Run Steps 1-4 in BilaterialRegistration-3D-Method1-PointMatching FIRST %%%%%%%%%%%%%%%%%%%%%

%% Display bilaterial path matches on the images

%pauses each match, so must press enter to see next match

mi1 = min(min(min(Im1)));
ma1 = max(max(max(Im1)));


mi2 = min(min(min(Im2)));
ma2 = max(max(max(Im2)));

for index1 = 1:numMatch
    tempIm1 = Im1;
    tempIm2 = Im2;
    
    pairedPointIm1 = optimalMatchFullDetailsBilat(index1,2:3);
    
    point1 = listPointIm1(pairedPointIm1(1),:); % point1 = [y1 x1]
    point2 = listPointIm1(pairedPointIm1(2),:); % point2 = [y2 x2]
    
    path1 = findPath3D(point1, point2, Im1, pathOpt,WoptPath,Vf1,2);

    
    pairedPointIm2 = optimalMatchFullDetailsBilat(index1,5:6);
    
    point1 = listPointIm2(pairedPointIm2(1),:);
    point2 = listPointIm2(pairedPointIm2(2),:);
       
    path2 = findPath3D(point1, point2, Im2, pathOpt,WoptPath,Vf2,2);

    %figure()
    figure('Name',strcat('Image1: Path ',int2str(index1)));
	
    x1 = path1(:,2);
    y1 = path1(:,1);
    z1 = path1(:,3);

    %%%%%%%%DISPLAY IMAGE AND PATH
    %isosurface(gm1,.05);axis equal;alpha(0.05)
    slice(gm1, [], [size(gm1,1)/2], [z1(1) z1(size(x1, 1))]) %, size(gm1,3)/2)
    %--------------------------old----------------isosurface(Im1,0);axis equal; alpha(0.5)  
     hold on;
     plot3(x1,y1,z1,'r*-')
     pause
    figure()
    
    x2 = path2(:,2);
    y2 = path2(:,1);
    z2 = path2(:,3);

    %%%%%%%%DISPLAY IMAGE AND PATH
    isosurface(gm2,.05);axis equal;alpha(0.05) 
    %isosurface(Im2,0);axis equal
    % alpha(0.5)  
     hold on;
     plot3(x2,y2,z2,'r-')%,'*-')
    
    pause;
end
