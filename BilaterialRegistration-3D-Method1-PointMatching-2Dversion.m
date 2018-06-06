%% Step 1: Set Path and Load your images

addpath('3D-HelperFunctions')
addpath('3D-HelperFunctions\FastMarching_version3b\')
addpath('3D-HelperFunctions\FastMarching_version3b\functions')
addpath('3D-HelperFunctions\FastMarching_version3b\shortestpath')
addpath('3D-HelperFunctions\toolbox_fast_marching')
addpath('3D-HelperFunctions\toolbox_fast_marching\toolbox')
addpath('3D-HelperFunctions\toolbox_fast_marching\tests')

Im1= phantom3d('modified shepp-logan');
Im1=Im1(:,:,25);

Im2 = imrotate(Im1, 30, 'bicubic', 'crop');


%% Step 2: Set variable that are needed later

%variables in Step 3
t1 = [0.1 0.5 0.1]; %threshold used when finding the feature points
t2 = [0.1 0.5 0.1]; %threshold used when finding the feature points
numFeat = 50; %the number of feature points when finding them with getDistantStrongestPoints2


%variables in Step 4
numBin = 4; %number of bins for the histogram feature descriptor (NOT used for 3D)
numEllipse=5; %number of feature descriptor for the ellipse around the path (NOT used for 3D)
pathOpt = 2; %path option in findPath3D
avgBlockPoint = 5; %if you want 5 or 3 box around your start and end point (must be 3 or 5)

%WoptPath=4;
WoptPath=2; %another option when using findPath3D

%numResample = 100; %number of points that you want to resample from the path(note it will be numResample)
    %numResample is calculated later on, it will be the number of points
    %in the largest path for all bilateral points
numBreaks = 8; % number of breaks you want in your path (how many sections you want in your path)

numMatch = 10; % the number of bilaterial path matches between the two images


%% Step 3: Option 2 - Find Feature Point (Find the gradient of the image)
%This option uses only gradient of the images to retrieve information

%gradient
[vx1 vy1] = gradient(Im1);
gm1=sqrt(vx1.^2+vy1.^2);
gm1 = mat2gray(gm1);

[vx2 vy2] = gradient(Im2);
gm2=sqrt(vx2.^2+vy2.^2);
gm2 = mat2gray(gm2);

S1=zeros(1,size(gm1,1),size(gm1,2));
fprintf('Start to find feature points in Image 1\n');

S1(1,:,:) = gm1;

[P1] = getDistantStrongestPoints2( S1 , numFeat , 0, t1);

fprintf('Completed finding the feature points in Image 1\n');

%find feature points for image 2

fprintf('Start to find feature points in Image 2\n');

S2=zeros(1,size(gm2,1),size(gm2,2));

S2(1,:,:) = gm2;

[P2] = getDistantStrongestPoints2(  S2 , numFeat , 0, t2);

fprintf('Completed finding the feature points in Image 2\n');

listPointIm1 = P1;
listPointIm2 = P2;


%% select a few points for testing

%listPointIm1 = [listPointIm1(30,:);listPointIm1(31,:);listPointIm1(35,:);listPointIm1(36,:);listPointIm1(37,:);listPointIm1(38,:)];
%listPointIm1(1,:)= [61,206,249]
%listPointIm1(5,:)= [17,104,95]
%listPointIm2 = [listPointIm2(30,:);listPointIm2(33,:);listPointIm2(34,:); 10, 219, 249; 3,136,136;31,100,95];


%% Step 4: Find the Feature Descriptors

%Feature descriptors that are used in 3D are
%   i) the euclidean distance along the path
%   ii)average pixel intensity for each section along the path (the number of sections depends on the number of breaks you choose)
%   iii) average pixel of a 3x3x3 or 5x5x5 around the start and end point of the path 
%   iv)curvature of the path

fprintf('Start to find feature descriptors in Image 1\n');


%find max length of paths between all points, so that we can resample all
%the paths to the largest path number (this is to minimize the problem of
%resampling points)

numResample1 = 0;

for index1 = 2:size(listPointIm1,1)
    index1
    point1 = listPointIm1(index1,:); %point1(2) is x value and point1(1) is y value and point1(3) is z value
    for index2 = 1:size(listPointIm1,1)
        if (index2 < index1)
            %remember: all the matrices are squared so row ith and jth column represents
            %starting point is i and the end point is j
            %   so there will be redundency so we want everything below the
            %   diagonal or above it
            point2 = listPointIm1(index2,:); %point2(2) is x value and point2(1) is y value and point1(3) is z valuea
           %find path
            path = findPath3D(point1, point2, Im1, pathOpt, WoptPath, gm1,2);
            size(path,1)
            if numResample1 < size(path,1)
               numResample1 = size(path,1);
            end
        end
    end
end

numResample2 = 0;

for index1 = 2:size(listPointIm2,1)
    index1
    point1 = listPointIm2(index1,:); %point1(2) is x value and point1(1) is y value and point1(3) is z value
    for index2 = 1:size(listPointIm2,1)
        if (index2 < index1)
            %remember: all the matrices are squared so row ith and jth column represents
            %starting point is i and the end point is j
            %   so there will be redundency so we want everything below the
            %   diagonal or above it
            point2 = listPointIm2(index2,:); %point2(2) is x value and point2(1) is y value and point1(3) is z valuea
           %find path
            path = findPath3D(point1, point2, Im2, pathOpt, WoptPath, gm2,2);
            size(path,1)
            if numResample2 < size(path,1)
               numResample2 = size(path,1);
            end
        end
    end
end

numResample = max(numResample1, numResample2);
%%
%create all the matrixes to store the feature descriptor values once they
%are found for Image 1

%NOTE: all the matrices are squared so row ith and jth column represents
%starting point is i and the end point is j

avgIPathsIm1 = zeros(size(listPointIm1,1), size(listPointIm1,1)); %average intensity along the whole path, not used since it is not helpful
lenPathsIm1 = zeros(size(listPointIm1,1), size(listPointIm1,1));
histPathsIm1 = zeros(size(listPointIm1,1), size(listPointIm1,1), numBin);
ellipsePathsIm1 = zeros(size(listPointIm1,1), size(listPointIm1,1),numEllipse);
avgPoint1Im1 = zeros(size(listPointIm1,1), size(listPointIm1,1));
avgPoint2Im1 = zeros(size(listPointIm1,1), size(listPointIm1,1));
avgIBreakPathsIm1 = zeros(size(listPointIm1,1), size(listPointIm1,1),numBreaks);

%find out the number of bilaterial paths that we will be working with
 numPaths=0;
 for index1 = 2:size(listPointIm1,1)
     for index2 = 1:size(listPointIm1,1)
         if (index2 < index1)
            numPaths = numPaths+1; 
         end
     end
 end
 
 %keep track of paths between points that are not successful (for some reason 
 %the fast marching sometimes creates paths with only one points and I can not 
 %figure out why so I was documenting when this happened)
 
 failedPathIm1 = 0;
 failedPathStartIm1 = [];
 failedPathEndIm1 = [];
 
allCurvIm1 = zeros(numResample-1, numPaths);

currPathNum = 1;
 
for index1 = 2:size(listPointIm1,1)
    index1
    point1 = listPointIm1(index1,:); %point1(2) is x value and point1(1) is y value and point1(3) is z value
    
    for index2 = 1:size(listPointIm1,1)

        
        if (index2 < index1)
            %remember: all the matrices are squared so row ith and jth column represents
            %starting point is i and the end point is j
            %   so there will be redundency so we want everything below the
            %   diagonal or above it
            
            point2 = listPointIm1(index2,:); %point2(2) is x value and point2(1) is y value and point1(3) is z valuea

            %calculate the avg pixel(3x3 or 5x5 box) of the start and end point  
            
            avgPoint1Im1(index1, index2) = pointAvg3D(Im1, point1, avgBlockPoint);
            avgPoint2Im1(index1, index2) = pointAvg3D(Im1, point2, avgBlockPoint);
            
           %find path
             
            path = findPath3D(point1, point2, Im1, pathOpt, WoptPath, gm1,2);
            
            %see if the findPath failed (ie. it says only one point is in the path which is not true)
            if size(path,1)==1
                 failedPathIm1 = failedPathIm1 + 1;
                 failedPathStartIm1 = [failedPathStartIm1; point1];
                 failedPathEndIm1 = [failedPathEndIm1; point2];
            else
                %resample the path
                x = path(:,2);
                y = path(:,1);
                z = path(:,3);
%                acc = accum(x,y,z);
%                resampPath = resampleLine3D2( x', y', z', acc, numResample);
                %fprintf( 'calc resampl')
                resampPath = resampleLine3Dsimpler( x, y, z, numResample);
                %resample is in hte order [y x z]
             %   plot3(x,y,z)
              %  hold on;
               % plot3(resampPath(:,2), resampPath(:,1), resampPath(:,3), 'r')
                 
                
                %calculate the curvature of the path using the new sampled
                %points that are evenly spaced out
                currCurv = curvM(resampPath(:,2),resampPath(:,1),resampPath(:,3));
                allCurvIm1(:, currPathNum) = currCurv';

                %calculate length of the pathway
                lenPathsIm1(index1, index2) = lengthPath3D(path);
                
                %calculate the average pixel intensity along the path for
                %each section according to the choosen breaks
                
                current = 1;
                T = findBreak(path, numBreaks);

                for bp = 1:numBreaks
                    if (T(bp)==0)
                        avgIBreakPathsIm1(index1, index2, bp) = 0;
                    else
                        avgIBreakPathsIm1(index1, index2, bp) = avgPathInt(path(current:current+(T(bp)-1),:),Im1);
                        current = current + T(bp);  
                    end  
                end


                currPathNum = currPathNum + 1;
            end
        end
    end
end

fprintf('Start to find feature descriptors in Image 2\n');

%create all the matrixes to store the feature descriptor values once they
%are found for Image 2

%NOTE: all the matrices are squared so row ith and jth column represents
%starting point is i and the end point is j

avgIPathsIm2 = zeros(size(listPointIm2,1), size(listPointIm2,1));%average intensity along the whole path, not used since it is not helpful
lenPathsIm2 = zeros(size(listPointIm2,1), size(listPointIm2,1));
histPathsIm2 = zeros(size(listPointIm2,1), size(listPointIm2,1), numBin);
ellipsePathsIm2 = zeros(size(listPointIm2,1), size(listPointIm2,1),numEllipse);
avgPoint1Im2 = zeros(size(listPointIm2,1), size(listPointIm2,1));
avgPoint2Im2 = zeros(size(listPointIm2,1), size(listPointIm2,1));
avgIBreakPathsIm2 = zeros(size(listPointIm2,1), size(listPointIm2,1),numBreaks);
invAvgIBreakPathsIm2 = zeros(size(listPointIm2,1), size(listPointIm2,1),numBreaks);

%find out the number of bilaterial paths that we will be working with
 numPaths=0;
 for index1 = 2:size(listPointIm2,1)
     for index2 = 1:size(listPointIm2,1)
         if (index2 < index1)
            numPaths = numPaths+1; 
         end
     end
 end
 
 %keep track of paths between points that are not successful (for some reason 
 %the fast marching sometimes creates paths with only one points and I can not 
 %figure out why so I was documenting when this happened)
 failedPathIm2 = 0;
 failedPathStartIm2 = [];
 failedPathEndIm2 = [];
 
 allCurvIm2 = zeros(numResample-1, numPaths);

 currPathNum = 1;
 
for index1 = 1:size(listPointIm2,1)
    index1
    point1 = listPointIm2(index1,:);
    
    for index2 = 1:size(listPointIm2,1)
        if (index2 < index1)
            
            %remember: all the matrices are squared so row ith and jth column represents
            %starting point is i and the end point is j
            %   so there will be redundency so we want everything below the
            %   diagonal or above it
            point2 = listPointIm2(index2,:); % point1 and point2 are in formation [y x]

            %calculate the avg pixel(3x3 or 5x5 box) of the start and end point 

            avgPoint1Im2(index1, index2) = pointAvg3D(Im2, point1, avgBlockPoint);
            avgPoint2Im2(index1, index2) = pointAvg3D(Im2, point2, avgBlockPoint);
                
            %find path
            path = findPath3D(point1, point2, Im2, pathOpt,WoptPath, gm2,2);
            
            if size(path,1)==1
                 failedPathIm2 = failedPathIm2 + 1;
                 failedPathStartIm2 = [failedPathStartIm2; point1];
                 failedPathEndIm2 = [failedPathEndIm2; point2];
      
            else
                %resample the points along the path
                
                x = path(:,2);
                y = path(:,1);
                z = path(:,3);
%                acc = accum(x,y,z);
%                resampPath = resampleLine3D( x', y', z', acc, numResample);
                resampPath = resampleLine3Dsimpler( x, y, z, numResample);

                %calculate the curvature of the path using the new sampled
                %points that are evenly spaced out
%                currCurv = curvM(resampPath(:,2),resampPath(:,1),resampPath(:,3));
%                allCurvIm2(:, currPathNum) = currCurv(1:numResample-2)';
                
                currCurv = curvM(resampPath(:,2),resampPath(:,1),resampPath(:,3));
                allCurvIm2(:, currPathNum) = currCurv';

                 %calculate length of the pathway
                lenPathsIm2(index1, index2) = lengthPath3D(path);

                %calculate the average pixel intensity along the path for
                %each section according to the choosen breaks
                current = 1;
                T = findBreak(path, numBreaks);
                for bp = 1:numBreaks
                    if (T(bp)==0)
                        avgIBreakPathsIm2(index1, index2, bp) = 0;
                    else
                        avgIBreakPathsIm2(index1, index2, bp) = avgPathInt(path(current:current+(T(bp)-1),:),Im2);
                        current = current + T(bp);  
                    end  
                end

                %calculate the average pixel intensity for the reverse path for
                %each section according to the choosen breaks
                    %this is need because the start point of image 1 my
                    %match with the end point of image 2 and the end point
                    %of image 1 may match with start point in image 2
                current = 1;
                T = T(end:-1:1);
                for bp = 1:numBreaks
                    if (T(bp)==0)
                        invAvgIBreakPathsIm2(index1, index2, bp) = 0;
                    else
                        invAvgIBreakPathsIm2(index1, index2, bp) = avgPathInt(path(current:current+(T(bp)-1),:),Im2);
                        current = current + T(bp);  
                    end  
                end
                
                currPathNum = currPathNum + 1;

            end
            
        end
        
    end
end




%Take all the matrices for the descriptors of each bilaterial paths and
%combine them into one matrix so we can manage them easier (for image 1)
%the final matrix is called allComboPointIm1

    %allCombo: the first column is start point of the path, column 2 is end
    %point for each column and column 3 is the average intensity at start
    %point, column 4 is the average intensity at end point and ...
    
    %now row 1 corresponds to bilaterial path 1 in image 1, row 2 is the
    %second bilaterial path in image 1, .... and the nth row is the nth
    %bilaterial path in image 1
    
allComboPointIm1 = [];

for index = 2:size(listPointIm1,1)

    point1Full = transpose(1:1:(index-1));
    point2Full = zeros(index - 1,1) + index;
    
    avgP = avgIPathsIm1(index, 1:index-1);
    avgL = lenPathsIm1(index, 1:index-1);
    his = histPathsIm1(index, 1:index-1, :);
    ell = ellipsePathsIm1(index, 1:index-1, :);
    avgP1 = avgPoint1Im1(index, 1:index-1);
    avgP2 = avgPoint2Im1(index, 1:index-1);
    avgBre = avgIBreakPathsIm1(index, 1:index-1, :);
    
    for i =1:size(point1Full,1)
        %need the second for loop because remember the matrices that are
        %store are the number of feature points in image 1 by number of
        %feature points in image one and there is a redundency since we
        %only want the values between the diagonal 
        point1 = point1Full(i);
        point2 = point2Full(i);
        
        addR = [point1, point2, avgP(i), avgL(i)];

        for b = 1:numBin
            addR = [addR, his(1,i,b)];
        end

        for e = 1:numEllipse
            addR = [addR, ell(1,i,e)];
        end

        addR = [addR, avgP1(i), avgP2(i)];
        addR = [addR, avgP2(i), avgP1(i)];

        for bp = 1:numBreaks
            addR = [addR, avgBre(1,i,bp)];
        end

        for bp = 1:numBreaks
            addR = [addR, avgBre(1,i,bp)];
        end

        allComboPointIm1 = [allComboPointIm1; addR];

    end
end

%Take all the matrices for the descriptors of each bilaterial paths and
%combine them into one matrix so we can manage them easier (for image 2).
%the final matrix is called allComboPointIm2

    %allCombo: the first column is start point of the path, column 2 is end
    %point for each column and column 3 is the average intensity at start
    %point, column 4 is the average intensity at end point and ...
    
    %now row 1 corresponds to bilaterial path 1 in image 2, row 2 is the
    %second bilaterial path in image 2, .... and the nth row is the nth
    %bilaterial path in image 2
    
allComboPointIm2 = [];
for index = 2:size(listPointIm2,1)
    point1Full = transpose(1:1:(index-1));
    point2Full = zeros(index - 1,1) + index;
    
    avgP = avgIPathsIm2(index, 1:index-1);
    avgL = lenPathsIm2(index, 1:index-1);
    his = histPathsIm2(index, 1:index-1, :);
    ell = ellipsePathsIm2(index, 1:index-1, :);
    avgP1 = avgPoint1Im2(index, 1:index-1);
    avgP2 = avgPoint2Im2(index, 1:index-1);
    avgBre = avgIBreakPathsIm2(index, 1:index-1, :);
    invAvgBre = invAvgIBreakPathsIm2(index, 1:index-1, :);
    
    for i =1:size(point1Full,1)
        %need the second for loop because remember the matrices that are
        %store are the number of feature points in image 2 by number of
        %feature points in image one and there is a redundency since we
        %only want the values between the diagonal 
        point1 = point1Full(i);
        point2 = point2Full(i);
        
        addR = [point1, point2, avgP(i), avgL(i)];

        for b = 1:numBin
            addR = [addR, his(1,i,b)];
        end

        for e = 1:numEllipse
            addR = [addR, ell(1,i,e)];
        end

        addR = [addR, avgP1(i), avgP2(i)];
        addR = [addR, avgP2(i), avgP1(i)];

        for bp = 1:numBreaks
            addR = [addR, avgBre(1,i,bp)];
        end

        for bp = 1:numBreaks
            addR = [addR, invAvgBre(1,i,bp)];
        end

        allComboPointIm2 = [allComboPointIm2; addR];
    end
end

FV1 = allComboPointIm1(:,3:size(allComboPointIm1,2));
FV2 = allComboPointIm2(:,3:size(allComboPointIm2,2));

nonNFV1 = FV1;
nonNFV2 = FV2;

%Normalize each feature

%normalize separately
% for j = 1:size(FV1,2)
%     FV1(:,j) = mat2gray(FV1(:,j));
%     FV2(:,j) = mat2gray(FV2(:,j));
% end



%normalize together
tempFVtog = [FV1; FV2];
for j = 1:size(FV1,2)
    tempFVtog(:,j) = mat2gray(tempFVtog(:,j));
end

FV1 = tempFVtog(1:size(FV1,1), 1:size(FV1,2));
FV2 = tempFVtog((size(FV1,1) + 1): (size(FV1,1) + size(FV2,1)), 1:size(FV2,2));

%% Step 4: Create the weights for each feature descriptor AND Compare bilaterial paths between images
%compare each bilaterial path in image 1 with all bilaterial path in image 2 (get a value describing how different they are )


%create weights for feature descriptors
ww = [];
ww = [ww 0];            %average intensity
ww = [ww 4];            %length
m = zeros(1, numBin) + 0;
ww = [ww m];      %four bins of histogram - DO NOT CHANGE THIS
ww = [ww 0];            %ellipse MinorAxisLength;
ww = [ww 0];            %ellipse MajorAxisLength;
ww = [ww 0];            %ellispe minor / major ;
ww = [ww 0];            %ellipse Eccentricity;
ww = [ww 0];            %ellipse Perimeter;
ww = [ww 0.6667];       %3x3 or 5x5 avg intensity start point
ww = [ww 0.6667];       %3x3 or 5x5 avg intensity end point
ww = [ww 0.6667];       %3x3 or 5x5 avg intensity start point
ww = [ww 0.6667];       %3x3 or 5x5 avg intensity end point (switch)
m = zeros(1, numBreaks) + 0.4;
ww = [ww m];            %average intensity broken into 8
ww = [ww m];            %average intensity broken into 8 (other way for Image 2)

%weight for curvature has to be separate
Wcurv = .3;


%create a matrix with all zeros except the diagonal values so identity
%matrix but the values on the diagonal will be our weight values
W = makeW(ww);

%Curvature comparison is a special case so we do this first and incorporate
%it later
curvCompare = zeros(size(allCurvIm1,2),size(allCurvIm2,2));
curvCompareOpp = zeros(size(allCurvIm1,2),size(allCurvIm2,2));

for mm = 1:size(allCurvIm1,2)
    %image 1 (rows)
    curv1 = allCurvIm1(:, mm);
    
    for nn = 1: size(allCurvIm2,2)
        %image 2 (cols)
        curv2 = allCurvIm2(:, nn);
        
        %in case they are opposites        
        curv2opp = curv2(size(curv2,1):-1:1);

        %sum over all the difference in the curvature and times it by the
        %weight
        curvCompare(mm,nn) = sum(abs(curv1-curv2))*Wcurv;
        curvCompareOpp(mm,nn) = sum(abs(curv1-curv2opp))*Wcurv;

    end    
end

FV1 = abs(FV1);
FV2 = abs(FV2);
diff = zeros(size(FV1,1), size(FV2,1));

%compare the rest of the descriptors with our previously found curvature
%difference
for mm = 1:size(FV1,1)
    %image 1 (rows)
    for nn = 1: size(FV2,1)
        %image 2 (cols)
        
        %compare the rest of the descriptors
        FVw = abs(FV1(mm,:) - FV2(nn,:))*W;
        
        %combine this difference with our curvature difference
        diff(mm,nn) = sum(FVw) + curvCompare(mm,nn) + curvCompareOpp(mm,nn); 
    end
    
end

%now we have to compensate for adding both assuming
%   1) start point in image 1 matches with start point in image 2, and end
%   point in image 1 matches end point in image 2
%   2) start point in image 1 matches with end point in image 2, and end
%   point in image 1 matches start point in image 2
%this bad so we must figure out which one give the lowest difference and
%then substract the difference that isnt being used

diff = abs(diff);

indDiff = zeros(size(diff));
diffDiff1 = zeros(size(diff));
diffDiff2 = zeros(size(diff));
%12,13,14,15
for jj = 1:size(diff, 1)
    %row jj (Image1)
    for kk = 1:size(diff, 2)
        %col kk (Image2)
        
        %start and end point in image 1 (use the average )
        S1 = FV1(jj,8+numBin);
        E1 = FV1(jj,9+numBin);
        %start and end point in image 2
        S2 = FV2(kk,8+numBin);
        E2 = FV2(kk,9+numBin);

        %now choose the split average intensity of the path
        endOfB1 = (12+numBin+numBreaks)-1;
        startOfB2 = (12+numBin+numBreaks);
        endOfB2 = (startOfB2+numBreaks)-1;

        diffDiff1(jj,kk) = (abs(S1-S2) + abs(E1-E2));  
        diffDiff2(jj,kk) = (abs(S1-E2) + abs(E1-S2));
        
        if ((abs(S1-S2) + abs(E1-E2)) < (abs(S1-E2) + abs(E1-S2)))
            %we want to keep the difference between images when start point
            %of image 1 matches closest to start point in image 2 and end
            %point in image 1 matches closest to end point in image 2 
            %AND
            %subtract the other difference between images when start point
            %of image 1 matches closest to end point in image 2 and end
            %point in image 1 matches closest to start point in image 2
            
            indDiff(jj,kk) = 1;
            
            breakCompensation = sum(abs((FV1(jj,startOfB2:endOfB2) - FV2(kk,startOfB2:endOfB2))*W(startOfB2:endOfB2,startOfB2:endOfB2)));
            startEndCompensation = sum(abs((FV1(jj,10+numBin:11+numBin) - FV2(kk,10+numBin:11+numBin))*W(10+numBin:11+numBin,10+numBin:11+numBin)));
            
            mustSub = breakCompensation + startEndCompensation;
            
             diff(jj,kk) = abs(diff(jj,kk) - mustSub);

             %subtract curv now since it is separate
             
             diff(jj,kk) = abs(diff(jj,kk) - curvCompareOpp(jj,kk));
             
        else
            %subtract the other difference between images when start point
            %of image 1 matches closest to end point in image 2 and end
            %point in image 1 matches closest to start point in image 2
            %AND           
            %we want to keep the difference between images when start point
            %of image 1 matches closest to start point in image 2 and end
            %point in image 1 matches closest to end point in image 2 
            
            indDiff(jj,kk) = 2;
            
            breakCompensation = sum(abs((FV1(jj,12+numBin:endOfB1) - FV2(kk,12+numBin:endOfB1))*W(12+numBin:endOfB1, 12+numBin:endOfB1)));
            startEndCompensation = sum(abs((FV1(jj,8+numBin:9+numBin) - FV2(kk,8+numBin:9+numBin))*W(8+numBin:9+numBin,8+numBin:9+numBin)));
            
            mustSub = breakCompensation + startEndCompensation;
            
            diff(jj,kk) = abs(diff(jj,kk) - mustSub);
            
            diff(jj,kk) = abs(diff(jj,kk) - curvCompare(jj,kk));
        end
    end
end

diff = abs(diff);

%not that we have the differences between all bilaterial paths between
%images, we can now match them
numMatch = 3
[optimalMatchBilat optimalMatchFullDetailsBilat] = matchTwoSetsBilat(diff, numMatch, allComboPointIm1, allComboPointIm2, 1);




%NOW use the other files to display various results




%% Save matrices in case you want to use them again

%save feature points
dlmwrite('storeListPointIm1-april26.txt', listPointIm1, 'delimiter', '\t', 'precision', 6);
dlmwrite('storelistPointIm2-april26.txt', listPointIm2, 'delimiter', '\t', 'precision', 6);

%save feature descriptors .... 
dlmwrite('storeFV1-3D-S1S2-2.txt', FV1, 'delimiter', '\t', 'precision', 6);
dlmwrite('storeFV2-3D-S1S2-2.txt', FV2, 'delimiter', '\t', 'precision', 6);

dlmwrite('allComboIm1-S1S2-2.txt', allComboPointIm1, 'delimiter', '\t', 'precision', 6);
dlmwrite('allComboIm2-S1S2-2.txt', allComboPointIm2, 'delimiter', '\t', 'precision', 6);

%% Load matrices

%load feature points
listPointIm1 = dlmread('storeListPointIm1-april26.txt');
listPointIm2 = dlmread('storeListPointIm2-april26.txt');

%load feature descriptors ... 
FV1 = dlmread('storeFV1-3D-S1S2.txt');
FV2 = dlmread('storeFV2-3D-S1S2.txt');

allComboPointIm1 = dlmread('allComboIm1-S1S2.txt');
allComboPointIm2 = dlmread('allComboIm2-S1S2.txt');