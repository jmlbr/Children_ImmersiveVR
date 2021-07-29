% Detect_discontinuities.m

% Script to indicate parts of kinematic data (head and torso angles) with
% discontinuities due to sensor jumps
% The waypoints containing these data points will not be considered
% for metrics based on kinematic data

% The codes identifies consecutive samples with angle differences larger
% than 20°. Pairs consisting of the index before and after the problematic
% index are written into a text file with the same name identifiers as the
% other data files. This file is then loaded by the script
% create_dataTable.m before computing the different kinematic metrics. 

% A prompt offers the option to enter the indices between which to reject the data 
% manually. This is useful in case several consecutive jumps occured

% ©Jenifer Miehlbradt, EPFL, 2021

%% 
dataFolder = cd; % ADD FOLDER CONTAINING DATA HERE;

manoeuvreListFiles = dir('s*manoeuvreList*');
manoeuvreListFiles = {manoeuvreListFiles.name};

droneFiles = dir('s*TimeAngleRatePosRot*');
droneFiles = {droneFiles.name};

WPdistFiles = dir('s*waypointDistTime*');
WPdistFiles = {WPdistFiles.name};

paramFiles = dir('s*miscExpPara*');
paramFiles = {paramFiles.name};

bodyAngleFiles = dir('s*BodyAngles*');
bodyAngleFiles = {bodyAngleFiles.name};

scoreFiles = dir('s*Score*');
scoreFiles = {scoreFiles.name};

%%

for ii = 1:length(droneFiles)
  
  % Find timestamp - part of the filename common to all file types from one
  % session
  nameParts = strsplit(droneFiles{ii},'_');
  timestamp = [nameParts{end-1},'_',nameParts{end}];
  currFiles = dir(['*',timestamp]);
  currFiles = {currFiles.name};
  
  % Load files
  droneVals = readDroneAngles(droneFiles{find(not(cellfun('isempty',strfind(droneFiles,timestamp))))});
  WPdist = importdata(WPdistFiles{find(not(cellfun('isempty',strfind(WPdistFiles,timestamp))))});
  bodyAngles = importdata(bodyAngleFiles{find(not(cellfun('isempty',strfind(bodyAngleFiles,timestamp))))});
  
  % Find when the drone started flying
  idx_start = find(diff(droneVals.data(:,end)) > 1,1,'first');
  idx_end = find(diff(droneVals.data(:,end)) < -1);
  if isempty(idx_end)
    idx_end = length(bodyAngles.data(:,3));
  end
  
  % Identify consecutive samples with angle differences > 20°
  if ~isempty(find(abs(diff(bodyAngles.data(1:idx_end,:)))> 20))...
      | ~isempty(find(abs(diff(droneVals.data(1:idx_end,2:4)))> 20))
    figure(1); clf
    subplot(2,1,1)
    plot(diff(droneVals.data(1:idx_end,2:4)))
    hold on
    y = ylim;
    plot([idx_start,idx_start], [y(1),y(2)],'k:')
    title('Difference drone angles')
    
    subplot(2,1,2)
    plot(diff(bodyAngles.data(1:idx_end,:)))
    hold on
    y = ylim;
    plot([idx_start,idx_start], [y(1),y(2)],'k:')
    title('Difference body angles')
    legend(bodyAngles.textdata)
    
    % Option to identify these indices automatically or to enter the values manually 
    answer = questdlg('Detect discontinuities automatically?', ...
      'Remove parts', ...
      'Yes','No','No');
    % Handle response
    switch answer
      case 'Yes'
        idx = [];
        for bodyAngleCol = 1:6
            idx = [idx;find(abs(diff(bodyAngles.data(1:idx_end,bodyAngleCol)))> 20)];
        end
        for droneAngleCol = 2:4
            idx = [idx;find(abs(diff(droneVals.data(1:idx_end,droneAngleCol)))> 20)];
        end
        idx = unique(idx);
        a = sort([idx-1;idx+1])';
        idx_to_remove = a;
      case 'No' % Enter pairs of values between which the data is to be excluded manually
        prompt = 'Idx to remove (start,stop): ';
        to_remove = input(prompt);
        idx_to_remove = to_remove;
    end
  else
    idx_to_remove = [];
  end
  
  % Write idx to remove to text file
  savenameParts = nameParts;
  savenameParts{5} = 'idxToRemove';
  saveName = strjoin(savenameParts,'_');
  dlmwrite(saveName,idx_to_remove)
  
  clear idx_to_remove
end





