clc
% cd ('/Users/miehlbra/Documents/Phd/Kids/Data/Recordings_Genova_Feb-March2018')
cd ('/Users/miehlbra/Documents/Phd/Kids/Data/Recordings_Genova')

files = dir('*.txt');
files = {files.name};

sub = {};
for ii = 1:9
  sub{end+1} = ['s00',num2str(ii)];
end
for ii = 10:24
  sub{end+1} = ['s0',num2str(ii)];
end

filesToCheck = {'_evaluationBefore_torso','_training1_torso','_training2_torso',...
  '_evaluationAfter_torso','_evaluationDayAfter_torso',...
  '_evaluationBefore_head','_training1_head','_training2_head',...
  '_evaluationAfter_head','_evaluationDayAfter_head'};

for s = [13,15:24]%length(sub)
  
  for ii = 1:length(filesToCheck)
    nFiles = length(find(not(cellfun('isempty',strfind(files,[sub{s},filesToCheck{ii}])))));
    if nFiles ~=6
      fprintf('Subject %s %s %d files \n',sub{s},filesToCheck{ii},nFiles)
      idx = find(not(cellfun('isempty',strfind(files,[sub{s},filesToCheck{ii}]))));
      %       for jj = 1:length(idx)
      %         disp(files{idx(jj)})
      %       end
      %       disp('..................................')
      
    end
  end
  
  
  
end

%% Delete unwanted files

clear all
close all
% clc

dataFolder = '/Users/miehlbra/Documents/Phd/Kids/Data/Recordings_Genova_Feb-March2018';
cd(dataFolder)
cd ('/Users/miehlbra/Documents/Phd/Kids/Data/Recordings_Genova')

manoeuvreListFiles = dir('s013*manoeuvreList*');
manoeuvreListFiles = {manoeuvreListFiles.name};

droneFiles = dir('s013*TimeAngleRatePosRot*');
droneFiles = {droneFiles.name};

WPdistFiles = dir('s013*waypointDistTime*');
WPdistFiles = {WPdistFiles.name};

paramFiles = dir('s013*miscExpPara*');
paramFiles = {paramFiles.name};

bodyAngleFiles = dir('s013*BodyAngles*');
bodyAngleFiles = {bodyAngleFiles.name};

scoreFiles = dir('s013*Score*');
scoreFiles = {scoreFiles.name};


% Delete unwanted files

for ii = 8:9%1:length(WPdistFiles)
  
  % Performance
  % Load distance to WP as computed by Unity
  WPdistTime = importdata(WPdistFiles{ii});
  
  % Load manoeuvre List file
  WPInfo = importdata(manoeuvreListFiles{ii});
  
  % Load corresponding TimeAngleRatePosRot file
  droneVals = readDroneAngles(droneFiles{ii});
  % find when the drone started flying
  idx_start = find(diff(droneVals.data(:,end)) > 1);
  trajCoord = droneVals.data(idx_start+1:end,8:10);
  RPYrates = droneVals.data(idx_start+1:end, 5:7);
  
  WPcoord = WPInfo.data;
  nWP = length(WPdistTime.data)-1;
  WPcoord = WPcoord(1:nWP,:);  % when unity decided to save 2 trajectories....
  
  % plot WP locations
  figure(20); clf; hold on
  plot3(WPInfo.data(1:nWP,1), WPInfo.data(1:nWP,3), WPInfo.data(1:nWP,2),'o','MarkerSize',10)
  %   axis equal
  xlabel('x')
  ylabel('z')
  zlabel('y')
  title(strrep(WPdistFiles{ii},'_',' '))
  view([-100, 33])
  
  % plot
  plot3(trajCoord(:,1),trajCoord(:,3),trajCoord(:,2))
  
  checkFile = questdlg('Keep file?','Check file','Yes','No', 'Yes');
  
  if strcmp (checkFile,'No')
    delete(manoeuvreListFiles{ii});
    delete(droneFiles{ii});
    delete(WPdistFiles{ii});
    delete(paramFiles{ii});
    delete(bodyAngleFiles{ii});
    delete(scoreFiles{ii});
    disp( ['deleting files ', manoeuvreListFiles{ii}])
  end
  
end



