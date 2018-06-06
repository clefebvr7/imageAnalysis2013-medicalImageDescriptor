%%%%%%%%%%%%%%%%%% Run Steps 1-4 in BilaterialRegistration-3D-Method1-PointMatching FIRST %%%%%%%%%%%%%%%%%%%%%

%% CODE 1: display all the feature points on the montage image and color code it according to the method uesd 

% also the montage image for each vesselness, gradient and phase symmetry is shown

addpath('3D-HelperFunctions')


montageRow = 27; %number of rows of images shown in the montage 
montageCol = 9; %number of columns of images shown in the montage 
%if you do not know then run CODE 2 below to observe it


rows = size(Im1,1);
cols = size(Im1,2);
N = floor(numFeat/size(S1,1)); %the number of points per image for feature selection

slices = size(Vf1,3);

coloursL=hsv(size(S1,1));
figure('Name','Vesselness');
montage(reshape(mat2gray(Vf1) , [rows cols 1 slices]))
title('Vesselness')
figure();
montage(reshape(mat2gray(psIm1) , [rows cols 1 slices]))
title('Phase symmetry')
figure();
montage(reshape(mat2gray(gm1) , [rows cols 1 slices]))
title('Gradient iamge')
figure();
montage(reshape(mat2gray(gm1) , [rows cols 1 slices]))
title('Original image with feature points')
hold on

colorInd = 1;
for i=1:size(listPointIm1,1)
    actRow = listPointIm1(i,1);
    actCol = listPointIm1(i,2);
    
    x = listPointIm1(i,2);
    y = listPointIm1(i,1);
    z = listPointIm1(i,3);
    
    while (z~=0)
        if z >= montageCol
            actRow = actRow + rows;
            z = z - montageCol;
        else 
            actCol = actCol + (z*cols);
            z = 0;
        end
    end
    
    if (i>colorInd*N)
        colorInd = colorInd + 1;
    end
    
    plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)

    pause
end

figure();
montage(reshape(mat2gray(gm2) , [rows cols 1 slices]))
title('Original image with feature points')
hold on

colorInd = 1;
%for i=1:size(listPointIm2,1)
for i=1:5
    actRow = listPointIm2(i,1);
    actCol = listPointIm2(i,2);
      
    x = listPointIm2(i,2);
    y = listPointIm2(i,1);
    z = listPointIm2(i,3);
    
    while (z~=0)
        if z >= montageCol
            actRow = actRow + rows;
            z = z - montageCol;
        else 
            actCol = actCol + (z*cols);
            z = 0;
        end
    end
    
    if (i>colorInd*N)
        colorInd = colorInd + 1;
    end
    
    plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)
    pause
end


%% CODE 1.2: display all the feature points on the montage image and color code it according to the method uesd 
%added different colors for specific points

% also the montage image for each vesselness, gradient and phase symmetry is shown

addpath('3D-HelperFunctions')


montageRow = 27; %number of rows of images shown in the montage 
montageCol = 9; %number of columns of images shown in the montage 
%if you do not know then run CODE 2 below to observe it


rows = size(Im1,1);
cols = size(Im1,2);
N = floor(numFeat/size(S1,1)); %the number of points per image for feature selection

slices = size(Vf1,3);


im1Match = optimalMatchBilat(:,1);
im2Match = optimalMatchBilat(:,2);

coloursL=hsv(size(im1Match,1)+1);

figure();
montage(reshape(mat2gray(gm1) , [rows cols 1 slices]))
title('Gradient image 1 with feature points')
hold on

colorInd = 1;
for i=1:size(listPointIm1,1)
    actRow = listPointIm1(i,1);
    actCol = listPointIm1(i,2);
    
    x = listPointIm1(i,2);
    y = listPointIm1(i,1);
    z = listPointIm1(i,3);
    
    while (z~=0)
        if z >= montageCol
            actRow = actRow + rows;
            z = z - montageCol;
        else 
            actCol = actCol + (z*cols);
            z = 0;
        end
    end
   
    plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)

end

counter = 0;

colorInd = colorInd + 1;

for pnt = 1:size(im1Match,1)
	for index1 = 2:size(listPointIm1,1)
		for index2 = 1:size(listPointIm1,1)
			if (index2 < index1)
				counter = counter + 1
				if counter == im1Match(pnt)	
					point1 = listPointIm1(index1,:);
					point2 = listPointIm1(index2,:);

					%do first point
					actRow = point1(1);
					actCol = point1(2);
    
					x = point1(2);
					y = point1(1);
					z = point1(3);
    				
    				while (z~=0)
    				    if z >= montageCol
    				        actRow = actRow + rows;
    				        z = z - montageCol;
    				    else 
    				        actCol = actCol + (z*cols);
        				    z = 0;
						end
    				end
	
					plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)
					
					%do second point
					actRow = point2(1);
					actCol = point2(2);
    
					x = point2(2);
					y = point2(1);
					z = point2(3);
    				
    				while (z~=0)
    				    if z >= montageCol
    				        actRow = actRow + rows;
    				        z = z - montageCol;
    				    else 
    				        actCol = actCol + (z*cols);
        				    z = 0;
						end
    				end
	
					plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)
					
					
					colorInd = colorInd + 1;
                end
			end
		end
	end
end

	
figure();
montage(reshape(mat2gray(gm2) , [rows cols 1 slices]))
title('Gradient image 2 with feature points')
hold on

colorInd = 1;
%for i=1:size(listPointIm2,1)
for i=1:5
    actRow = listPointIm2(i,1);
    actCol = listPointIm2(i,2);
    
    x = listPointIm2(i,2);
    y = listPointIm2(i,1);
    z = listPointIm2(i,3);
    
    while (z~=0)
        if z >= montageCol
            actRow = actRow + rows;
            z = z - montageCol;
        else 
            actCol = actCol + (z*cols);
            z = 0;
        end
    end
    
    plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)

end

counter = 0;

colorInd = colorInd + 1;

for pnt = 1:size(im2Match,1)
	for index1 = 2:size(listPointIm2,1)
		for index2 = 1:size(listPointIm2,1)
			if (index2 < index1)
				counter = counter + 1
				if counter == im2Match(pnt)	
					point1 = listPointIm2(index1,:);
					point2 = listPointIm2(index2,:);

					%do first point
					actRow = point1(1);
					actCol = point1(2);
    
					x = point1(2);
					y = point1(1);
					z = point1(3);
    				
    				while (z~=0)
    				    if z >= montageCol
    				        actRow = actRow + rows;
    				        z = z - montageCol;
    				    else 
    				        actCol = actCol + (z*cols);
        				    z = 0;
						end
    				end
	
					plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)
					
					%do second point
					actRow = point2(1);
					actCol = point2(2);
    
					x = point2(2);
					y = point2(1);
					z = point2(3);
    				
    				while (z~=0)
    				    if z >= montageCol
    				        actRow = actRow + rows;
    				        z = z - montageCol;
    				    else 
    				        actCol = actCol + (z*cols);
        				    z = 0;
						end
    				end
	
					plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)
					
					
					colorInd = colorInd + 1;
                end 
            end     
        end
    end
end

%% CODE 2:
%can't tell the number of images per row and column, so alternate images
%with black and white so u can see it better

test = zeros(size(Im1));

t = 1;
for i = 1:size(Im1,3)
    temp = zeros(size(Im1,1),size(Im1,2));
    
    if t==1
        temp = temp + 1;
        t=0;
    else
        
        t=1;
        
    end
    test(:,:,i) = temp;
end

        
figure();
montage(reshape(test , [rows cols 1 slices]))        
   





%% CODE 3: Display two points on the montage image
%for example given two points (start point s, end point e)

addpath('3D-HelperFunctions')

%CHANGE the s, e depending on which points you want to observe
s=[65 256 255];
e=[65 256 116];

montageRow = 27;
montageCol = 9;

rows = size(Im1,1);
cols = size(Im1,2);

slices = size(Vf1,3);

tempVf1 = thresholdIm ( Vf1, t, 0);

figure();
montage(reshape(mat2gray(tempVf1) , [rows cols 1 slices]))
hold on
coloursL=hsv(1);
SE = [s;e];
colorInd=1;
for i=1:size(SE,1)
    actRow = SE(i,1);
    actCol = SE(i,2);
    
    x = SE(i,2);
    y = SE(i,1);
    z = SE(i,3);
    
    while (z~=0)
        if z >= montageCol
            actRow = actRow + rows;
            z = z - montageCol;
        else 
            actCol = actCol + (z*cols);
            z = 0;
        end
    end

    
    plot(actCol, actRow, 'b.', 'color',[ coloursL(colorInd,1) coloursL(colorInd,2) coloursL(colorInd,3)], 'MarkerSize',20)

end
