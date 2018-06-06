%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Very Rough - Was in the process of converting from
% 2D to 3D and deciding to stick with method 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[dist, path, pred]=graphshortestpath(G,S,T)
addpath('3D-HelperFunctions')
%addpath('3D-HelperFunctions\sift')
%addpath('3D-HelperFunctions\livewire')

%load image
Im1= phantom3d('modified shepp-logan');
Im1 = Im1 * 100;

Im2 = imrotate(Im1, 30, 'bicubic', 'crop');

Im1=imfilter(Im1, fspecial('Gaussian',10,5), 'same');
Im2=imfilter(Im2, fspecial('Gaussian',10,5), 'same');


numBin=4;
numBreaks = 8;
numEllipse=5;

%% HARD CODED FEATURE POINT SELECTION

%listPointI2 = [listPointI2(1:2,:); 179 166; listPointI2(3,:); 179 166];

%listPointIm1 = [87 63 64; 52 55 52; 27 91 84; 64 64 115; 93 57 22; 64 64 15; 110 110 52]; %y, x, then z
listPointI2 = [86 77 64; 57 50 52; 22 75 84; 64 64 115; 92 73 22; 64 64 15; 100 20 84]; %y, x, then z

%% Compairson Matrix

numFeatures = 8;

angL = [0:10:360];
errors = zeros(size(angL));

%for Im2 determining the features for the paths between feature points is
%not in the for loop since we only do it once and image 2 never changes

avgIPathsI2 = zeros(size(listPointI2,1), size(listPointI2,1));
lenPathsI2 = zeros(size(listPointI2,1), size(listPointI2,1));
histPathsI2 = zeros(size(listPointI2,1), size(listPointI2,1), numBin);

avgPoint1I2 = zeros(size(listPointI2,1), size(listPointI2,1));
avgPoint2I2 = zeros(size(listPointI2,1), size(listPointI2,1));

ellipsePathsI2 = zeros(size(listPointI2,1), size(listPointI2,1),numEllipse);
avgIBreakPathsI2 = zeros(size(listPointI2,1), size(listPointI2,1),numBreaks);
invAvgIBreakPathsI2 = zeros(size(listPointI2,1), size(listPointI2,1),numBreaks);

numFeat = 1 + 1 + numBin + ((1 + 1)*2) + numEllipse + (numBreaks*2);

for index1 = 1:size(listPointI2,1)

    point1 = listPointI2(index1,:);

    for index2 = 1:size(listPointI2,1)

        if (index2 < index1)
            point2 = listPointI2(index2,:);
        
            %avg (3x3 or 5x5 box) end points
            avgPoint1I2(index1, index2) = pointAvg3D(I2, point1, 5);
            avgPoint2I2(index1, index2) = pointAvg3D(I2, point2, 5);
%------------------------------------------------------------------CONVERTED TO 3D-------------------------------------------------------------------------------------------

            
           %find path
            path = findPath(point1, point2, I2, 1,2);

            
            %calculate average pixel along the pathway
            avgIPath = characterForArea(path, I2);
            avgIPathsI2(index1, index2) = avgIPath;
            
            %calculate length along the pathway
            lenPathsI2(index1, index2) = lengthPath(path);
            
            %calculate histogram along the pathway
            
            bins = histPath(path, I2, numBin);
            
            for b = 1:numBin
                histPathsI2(index1, index2, b) = bins(b);
            end
            
            ellips = regProp(I2, path, numEllipse);
            
            for e = 1:numEllipse
                ellipsePathsI2(index1, index2, e) = ellips(e);
            end
            
            current = 1;
            
            T = findBreak(path, numBreaks);
            for bp = 1:numBreaks
                if (T(bp)==0)
                    avgIBreakPathsI2(index1, index2, bp) = 0;
                else
                    avgIBreakPathsI2(index1, index2, bp) = characterForArea(path(current:current+(T(bp)-1),:),I2);
                    current = current + T(bp);  
                end  
            end
            
            current = 1;
            T = T(end:-1:1);
            for bp = 1:numBreaks
                if (T(bp)==0)
                    invAvgIBreakPathsI2(index1, index2, bp) = 0;
                else
                    invAvgIBreakPathsI2(index1, index2, bp) = characterForArea(path(current:current+(T(bp)-1),:),I2);
                    current = current + T(bp);  
                end  
            end
            
        end
    end
end

maxS = (size(listPointI2,1)^2)/2;
storeFV_2 = zeros(maxS, numFeat, size(angL,2));
storeFVHypoIm = zeros(maxS, numFeat, size(angL,2));
indFV_2 = [];
indFVHypoIm = [];

for k = 1:size(angL,2)
    k
    
    I1 = phantom('Modified Shepp-Logan',200);
    hypoIm = imrotate(I1, angL(k), 'bicubic', 'crop');
    hypoIm=imfilter(hypoIm, fspecial('Gaussian',10,5), 'same');

    
    mask = imrotate(zeros(size(I1))+1, angL(k), 'bicubic', 'crop');

    filtFP = [];
    FPRemovedIndex = [];
    
  
    
    for m = 1:size(listPointI2,1)
        if (mask(listPointI2(m,1),listPointI2(m,2))>=0.5)
            %add feature points since 
            filtFP = [filtFP;listPointI2(m,:)];
        else
            FPRemovedIndex = [FPRemovedIndex; m];
        end
    end

   %remove paths that contain the removed feature points (have to do two loops to go through cols and then go through rows)
   tempFiltAvgIPathsI2 = [];
   tempFiltLenPathsI2 = [];
   tempFiltHistPathsI2 = [];
   
   tempAvgPoint1I2 = [];
   tempAvgPoint2I2 = [];
   
   tempFiltEllipsePathsI2 = [];
   tempAvgIBreakPathsI2 = [];
   tempInvAvgIBreakPathsI2 = [];
   
    for col = 1:size(listPointI2,1)
        if (~isMem(col,FPRemovedIndex))
            tempFiltAvgIPathsI2 = [tempFiltAvgIPathsI2 avgIPathsI2(:,col)];
            tempFiltLenPathsI2 = [tempFiltLenPathsI2 lenPathsI2(:,col)];
            tempFiltHistPathsI2 = [tempFiltHistPathsI2 histPathsI2(:,col,:)];

            tempAvgPoint1I2 = [tempAvgPoint1I2 avgPoint1I2(:,col,:)];
            tempAvgPoint2I2 = [tempAvgPoint2I2 avgPoint2I2(:,col,:)];
            
            tempFiltEllipsePathsI2 = [tempFiltEllipsePathsI2 ellipsePathsI2(:,col,:)];
            tempAvgIBreakPathsI2 = [tempAvgIBreakPathsI2 avgIBreakPathsI2(:,col,:)];
            tempInvAvgIBreakPathsI2 = [tempInvAvgIBreakPathsI2 invAvgIBreakPathsI2(:,col,:)];
        end
    end
    
   filtAvgIPathsI2 = [];
   filtLenPathsI2 = [];
   filtHistPathsI2 = [];
   
   filtAvgPoint1I2 = [];
   filtAvgPoint2I2 = [];
   
   filtEllipsePathsI2 = [];
   filtAvgIBreakPathsI2 = [];
   filtInvAvgIBreakPathsI2 = [];

    for row = 1:size(listPointI2,1)
         if (~isMem(row,FPRemovedIndex))
             filtAvgIPathsI2 = [filtAvgIPathsI2; tempFiltAvgIPathsI2(row,:)];
             filtLenPathsI2 = [filtLenPathsI2; tempFiltLenPathsI2(row,:)];
             filtHistPathsI2 = [filtHistPathsI2; tempFiltHistPathsI2(row,:,:)];
             
             filtAvgPoint1I2 = [filtAvgPoint1I2; tempAvgPoint1I2(row,:,:)];
             filtAvgPoint2I2 = [filtAvgPoint2I2; tempAvgPoint2I2(row,:,:)];
             
             filtEllipsePathsI2 = [filtEllipsePathsI2; tempFiltEllipsePathsI2(row,:,:)];
             filtAvgIBreakPathsI2 = [filtAvgIBreakPathsI2; tempAvgIBreakPathsI2(row,:,:)];
              
             filtInvAvgIBreakPathsI2 = [filtInvAvgIBreakPathsI2; tempInvAvgIBreakPathsI2(row,:,:)];

        end
    end
    
    %calculate their path descriptors for T(Image 1) and use filtFP for the
    %list of features points

    avgIPathsHypoIm = zeros(size(filtFP,1), size(filtFP,1));
    lenPathsHypoIm = zeros(size(filtFP,1), size(filtFP,1));
    histPathsHypoIm = zeros(size(filtFP,1), size(filtFP,1), numBin);
    
    avgPoint1HypoIm = zeros(size(listPointI2,1), size(listPointI2,1));
    avgPoint2HypoIm = zeros(size(listPointI2,1), size(listPointI2,1));

    ellipsePathsHypoIm = zeros(size(filtFP,1), size(filtFP,1),numEllipse);
    avgIBreakPathsHypoIm = zeros(size(filtFP,1), size(filtFP,1),numBreaks);

        
    for index1 = 2:size(filtFP,1)

        point1 = filtFP(index1,:);

        for index2 = 1:size(filtFP,1)

            if (index2 < index1)
                point2 = filtFP(index2,:);

                            
                %avg (3x3 or 5x5 box) end points
                avgPoint1HypoIm(index1, index2) = pointAvg(hypoIm, point1, 3);
                avgPoint2HypoIm(index1, index2) = pointAvg(hypoIm, point2, 3);
            
               %find path
                path = findPath(point1, point2, hypoIm, 1,2);

                %calculate average pixel along the pathway
                avgIPathsHypoIm(index1, index2) = characterForArea(path, hypoIm);

                %calculate length along the pathway
                lenPathsHypoIm(index1, index2) = lengthPath(path);

                %calculate histogram along the pathway

                bins = histPath(path, hypoIm, numBin);

                for b = 1:numBin
                    histPathsHypoIm(index1, index2, b) = bins(b);
                end
            
            
                ellips = regProp(hypoIm, path, numEllipse);

                for e = 1:numEllipse
                    ellipsePathsHypoIm(index1, index2, e) = ellips(e);
                end
            
                current = 1;

                            
                T = findBreak(path, numBreaks);
                for bp = 1:numBreaks
                    if (T(bp)==0)
                        avgIBreakPathsHypoIm(index1, index2, bp) = 0;
                    else
                        avgIBreakPathsHypoIm(index1, index2, bp) = characterForArea(path(current:current+(T(bp)-1),:),hypoIm);
                        current = current + T(bp);  
                    end  
                end

            end
        end
    end

    %isolate the actual points from the above paired point matrices
    %allCombo: the first column is value, column 2 is labeling point 1 and column 3 is point 2 where column2,3 are the combined points    
    allComboPointHypoIm = [];
    
    for index = 2:size(filtFP,1)
        point1 = transpose(1:1:(index-1));
        point2 = zeros(index - 1,1) + index;
        
        addR = [point1, point2, transpose(avgIPathsHypoIm(index, 1:index-1)), transpose(lenPathsHypoIm(index, 1:index-1))];

        for b = 1:numBin
            addR = [addR, transpose(histPathsHypoIm(index, 1:index-1, b))];
        end

        for e = 1:numEllipse
            addR = [addR, transpose(ellipsePathsHypoIm(index, 1:index-1, e))];
        end
    
        addR = [addR, transpose(avgPoint1HypoIm(index, 1:index-1)), transpose(avgPoint2HypoIm(index, 1:index-1))];
        addR = [addR, transpose(avgPoint2HypoIm(index, 1:index-1)), transpose(avgPoint1HypoIm(index, 1:index-1))];

        for bp = 1:numBreaks
            addR = [addR, transpose(avgIBreakPathsHypoIm(index, 1:index-1, bp))];
        end

        for bp = 1:numBreaks
            addR = [addR, transpose(avgIBreakPathsHypoIm(index, 1:index-1, bp))];
        end
        
        allComboPointHypoIm = [allComboPointHypoIm; addR];
    end

   

    allComboPointI_2 = [];
    for index = 2:size(filtFP,1)
        point1 = transpose(1:1:(index-1));
        point2 = zeros(index - 1,1) + index;
        
        addR = [point1, point2, transpose(filtAvgIPathsI2(index, 1:index-1)), transpose(filtLenPathsI2(index, 1:index-1))];

        for b = 1:numBin
            addR = [addR, transpose(filtHistPathsI2(index, 1:index-1, b))];
        end

        for e = 1:numEllipse
            addR = [addR, transpose(filtEllipsePathsI2(index, 1:index-1, e))];
        end 

        addR = [addR, transpose(avgPoint1I2(index, 1:index-1)), transpose(avgPoint2I2(index, 1:index-1))];
        addR = [addR, transpose(avgPoint2I2(index, 1:index-1)), transpose(avgPoint1I2(index, 1:index-1))];

        for bp = 1:numBreaks
            addR = [addR, transpose(filtAvgIBreakPathsI2(index, 1:index-1, bp))];
        end

        for bp = 1:numBreaks
            addR = [addR, transpose(filtInvAvgIBreakPathsI2(index, 1:index-1, bp))];
        end

        allComboPointI_2 = [allComboPointI_2; addR];
    end
    
   

    if (numFeatures == 1)
        FVHypoIm = allComboPointHypoIm(:,3);
        FV_2 = allComboPointI_2(:,3);
 
        %now weight the features and sum them up to get one value
        W = eye(1);
        
    else
         FVHypoIm = allComboPointHypoIm(:,3:size(allComboPointHypoIm,2));
         FV_2 = allComboPointI_2(:,3:size(allComboPointI_2,2));
         
         %now weight the features and sum them up to get one value
         W = eye(size(FV_2,2));
    end


    %normalize each feature
    for j = 1:size(FV_2,2)
        FVHypoIm(:,j) = mat2gray(FVHypoIm(:,j));
        FV_2(:,j) = mat2gray(FV_2(:,j));
    end
    
    %store FVs so dont have to keep doing it if just wants to changes the
    %W's
    storeFV_2(1:size(FV_2,1),:,k) = FV_2;
    storeFVHypoIm(1:size(FVHypoIm,1),:,k) = FVHypoIm;
    indFV_2 = [indFV_2 size(FV_2,1)];
    indFVHypoIm = [indFVHypoIm size(FVHypoIm,1)];
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    FVHypoImW = (FVHypoIm * W);
    DescrHypoIm = sum(FVHypoImW,2);

    FV_2W = (FV_2 * W);
    Descr_2 = sum(FV_2W,2);

    errors(k) = sum(sqrt((DescrHypoIm-Descr_2).^2));
end
    
    
figure(); plot(angL, errors);
title('Comparison between Angle of Transformation and the Total Error for 8 Features');
xlabel('Angle of Proposed Transformation');
ylabel('Total Error of Descriptors Between Images');
















