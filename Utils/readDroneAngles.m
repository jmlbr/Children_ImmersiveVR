function A = readDroneAngles(filename)

A = importdata(filename);
[trigIdx,~] = find(A.data(:,1) == 999);
%A.data = A.data(trigIdx(1):end,:);
[a,~] = find(A.data(:,1) == 999);
A.data(unique(a),:) = [];



% DroneAngles = [C{2},C{3},C{4}];
% % remove 999 entries
% [a,~] = find(DroneAngles == 999);
% trigIdx = unique(a);
% DroneAngles = DroneAngles(trigIdx(1):end,:);
% [a,] = find(DroneAngles == 999);
% DroneAngles(unique(a),:) = [];
% DroneAngles = resample(double(DroneAngles),100,250);

