%% display using slice

addpath('3D-HelperFunctions')


slice(mat2gray(Im1),150,40,150)

%% display using isosurface

figure
isosurface(mat2gray(Im1),0)
alpha(0.5)
hold on;
%plot3([1 50],[1 50],[1 50])

