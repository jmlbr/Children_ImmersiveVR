% createDataTable_flyExp.m

% Processing routine to extract kinematic parameters from the files
% recorded during the flight experiment.
% Variables:
% - distance to waypoint,
% - time between waypoints,
% - Head-torso correlations in all planes,
% - Head and torso movement amplitudes(quartiles),
% - head, torso and drone spectral arc length (SAL, smoothness metric),
% - head anchoring index in all planes,
% - peak of the head-torso cross correlation,
% - number of peaks in the velocity profile for the head, torso and drone movements,
% - mean and maximum torso and head velocities,
% - vertical and horizontal error (distance to waypoint),
% - dynamic time warp distance(DTW) between head and torso (similarity metric),
% - path ratio between optimal path and actual path between consecutive waypoints.

% COMPUTES ALL THE VARIABLES FOR SEGMENTS BETWEEN WAYPOINTS!

% Saves 2 data tables: one containing the values for each segment (T), one
% containing average values for each session (T_means)

% ©Jenifer Miehlbradt, EPFL, 2021

%% Modify paths here
clear all
close all
clc

% Helper functions
addpath('Utils')   % ADD 'Utils' FOLDER HERE;


mainFolder = cd % ADD FOLDER CONTAINING DATA HERE;
dataFolder = [mainFolder,'/rawFiles']; 
saveFolder = cd; % ADD FOLDER TO SAVE DATA TABLES HERE;
cd(dataFolder)

%% Get file lists
manoeuvreListFiles = dir('s*manoeuvreList*');
manoeuvreListFiles = {manoeuvreListFiles.name};

droneFiles = dir('s*TimeAngleRatePosRot*');
droneFiles = {droneFiles.name};

WPdistFiles = dir('s*waypointDistTime*');
WPdistFiles = {WPdistFiles.name};

bodyAngleFiles = dir('s*BodyAngles*');
bodyAngleFiles = {bodyAngleFiles.name};

idxToRemoveFiles = dir('s*idxToRemove*');
idxToRemoveFiles = {idxToRemoveFiles.name};


%% Load subject information
load([mainFolder,'/subjectInfo.mat']);
ages = subjectInfo.Age;
genders = subjectInfo.Gender;

% Adjustment to cluster age groups, remove if unnecessary
ages(find(ages == 11)) = 10; % to cluster 10 and 11 year olds
ages(find(ages > 20)) = 25;  % to cluster adults

%% Get unique subject identifiers
fileNameParts = cellfun(@(x) strsplit(x, '_'), manoeuvreListFiles, 'UniformOutput', false);
fileNameParts = vertcat(fileNameParts{:});
[~,idx]= (unique(cell2mat(fileNameParts(:,1)),'rows'));
subjects = fileNameParts(idx,1);
clear idx fileNameParts

%% Compute metrics and store data in table

% Initialize some variables
[ID,Age,Gender, Control, Phase, dist2WP,time,...
    HTcorr_pp, HTcorr_rr,HTcorr_ry,HTcorr_yr,HTcorr_yy,...
    HeadQuart_p, HeadQuart_r, HeadQuart_y,TorsoQuart_r, TorsoQuart_p, TorsoQuart_y,...
    torsoSAL, headSAL, droneSAL, AI_roll, AI_pitch, AI_yaw,...
    CCpeak_pp, CCpeak_rr,CCpeak_yy,CCpeak_ry,CCpeak_yr,HeadAngles, BodyAngles,...
    headPks_p, headPks_r, headPks_y, torsoPks_r,torsoPks_p, torsoPks_y,...
    dronePks_r,dronePks_p, dronePks_y,torsoCorr_rp, torsoCorr_ry, torsoCorr_py,...
    meanSpeed_torso_r,meanSpeed_torso_p,meanSpeed_torso_y,...
    maxSpeed_torso_r,maxSpeed_torso_p,maxSpeed_torso_y,...
    torsoSpeed_norm_mean, torsoSpeed_norm_max,...
    error_h, error_v,dtw_roll, dtw_pitch, dtw_yaw,path_ratio] = deal([]);

% Table to store the kinematic variables
T = table(ID,Age,Gender, Control, Phase, dist2WP,time,...
    HTcorr_pp, HTcorr_rr,HTcorr_ry,HTcorr_yr,HTcorr_yy,...
    HeadQuart_p, HeadQuart_r, HeadQuart_y,TorsoQuart_r, TorsoQuart_p, TorsoQuart_y,...
    torsoSAL, headSAL, droneSAL,AI_roll, AI_pitch, AI_yaw,...
    CCpeak_pp, CCpeak_rr,CCpeak_yy,CCpeak_ry,CCpeak_yr,...
    headPks_p, headPks_r, headPks_y, torsoPks_r,torsoPks_p, torsoPks_y,...
    dronePks_r,dronePks_p, dronePks_y,...
    torsoCorr_rp, torsoCorr_ry, torsoCorr_py,...
    meanSpeed_torso_r,meanSpeed_torso_p,meanSpeed_torso_y,...
    maxSpeed_torso_r,maxSpeed_torso_p,maxSpeed_torso_y,torsoSpeed_norm_mean, torsoSpeed_norm_max,...
    error_h, error_v,dtw_roll, dtw_pitch, dtw_yaw,path_ratio,...
    'VariableNames',{'ID','Age','Gender', 'Control', 'Phase', 'dist2WP','Time','HTcorr_pp', 'HTcorr_rr','HTcorr_ry','HTcorr_yr','HTcorr_yy',...
    'HeadQuart_p', 'HeadQuart_r', 'HeadQuart_y','TorsoQuart_r', 'TorsoQuart_p', 'TorsoQuart_y', 'TorsoSAL', 'HeadSAL', 'DroneSAL',...
    'AI_roll', 'AI_pitch', 'AI_yaw','CCpeak_pp', 'CCpeak_rr','CCpeak_yy','CCpeak_ry','CCpeak_yr',...
    'HeadPks_p', 'HeadPks_r', 'HeadPks_y', 'TorsoPks_r','TorsoPks_p', 'Torso_pks_y',...
    'DronePks_r','DronePks_p', 'DronePks_y','TorsoCorr_rp', 'TorsoCorr_ry', 'TorsoCorr_py',...
    'MeanSpeed_torso_r','MeanSpeed_torso_p','MeanSpeed_torso_y',...
    'MaxSpeed_torso_r','MaxSpeed_torso_p','MaxSpeed_torso_y','TorsoSpeed_norm_mean','TorsoSpeed_norm_max',...
    'Error_h', 'Error_v','Dtw_roll', 'Dtw_pitch', 'Dtw_yaw','Path_ratio'});

phases = {'Before','training1','training2','After', 'DayAfter'};
control_mode = ['h','t'];   % Head or Torso

SR = 1000/68;    % 1 sample all 68 ms =14.7 Hz

for ii = 1:length(subjects)
    
    for kk = 1:2
        
        for p = 1:5
            
            % Filenames
            if ismember(p, [1,4,5])
                common_name = [subjects{ii},'_evaluation',phases{p},'_',control_mode(kk)];
            elseif ismember(p, [2,3])
                common_name = [subjects{ii},'_',phases{p},'_',control_mode(kk)];
            end
            
            % Load files
            if ~isempty(find(not(cellfun('isempty',strfind(WPdistFiles,common_name)))))
                WPdistTime = importdata(WPdistFiles{find(not(cellfun('isempty',strfind(WPdistFiles,common_name))))});
                droneVals = readDroneAngles(droneFiles{find(not(cellfun('isempty',strfind(droneFiles,common_name))))});
                droneData = importdata(droneFiles{find(not(cellfun('isempty',strfind(droneFiles,common_name))))}); % needed to get the timestamps of WP crossing
                WPInfo = importdata(manoeuvreListFiles{find(not(cellfun('isempty',strfind(droneFiles,common_name))))}); %location of WP
                idxToRemove = importdata(idxToRemoveFiles{find(not(cellfun('isempty',strfind(idxToRemoveFiles,common_name))))}); % parts of the file to remove (drone + kin data)
                
                bodyAngles = importdata(bodyAngleFiles{find(not(cellfun('isempty',strfind(bodyAngleFiles,common_name))))});
                
                switch p
                    case 1      % Before
                        nWP = 26;
                    case 2      % Training 1
                        nWP = 42;
                    case 3      % Training 2
                        nWP = 42;
                    case 4      % After
                        nWP = 18;
                    case 5      % Day After
                        nWP = 26;
                end
                
                
                %%%% Error (distance to waypoint) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if( WPdistTime.data(end, end) == -1 | isnan(WPdistTime.data(end, end)) )
                    dist2WP = WPdistTime.data(1:end-1,1);
                elseif ( WPdistTime.data(end, end) ==  WPdistTime.data(end-1, end)+1)
                    dist2WP = WPdistTime.data(1:end,1);
                else
                    dist2WP = [];
                    disp('??')
                end
                
                
                if length(dist2WP) < nWP
                    dist2WP(end+1:nWP) = nan;
                else
                    dist2WP= dist2WP(1:nWP);
                end
                
                %%% Error direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                errorDir = abs(getErrorDirection(WPInfo, droneVals, dist2WP, nWP));
                if length(errorDir) < nWP
                    errorDir(end+1:nWP,:) = nan;
                else
                    errorDir= errorDir(1:nWP,:);
                end
                
                %%% Time from previous WP %%%%%%%%%%%%%%%%%
                if( WPdistTime.data(end, end) == -1 | isnan(WPdistTime.data(end, end)) )
                    time = WPdistTime.data(1:end-1,2);
                elseif ( WPdistTime.data(end, end) ==  WPdistTime.data(end-1, end)+1)
                    time = WPdistTime.data(1:end,2);
                else
                    time = [];
                    disp('??')
                end
                
                
                if length(time) < nWP
                    time(end+1:nWP) = nan;
                else
                    time= time(1:nWP);
                end
                time(1) = nan;      % First duration takes into account pause before start
                
                
                % Find when the drone started flying
                idx_start = find(diff(droneVals.data(:,end)) > 1,1,'first');
                idx_end = find(diff(droneVals.data(:,end)) < -1);
                if isempty(idx_end)
                    idx_end = length(bodyAngles.data(:,3));
                end

                % Find index when waypoints were crossed, corresponds to '999'
                % entries in the droneData matrix
                idx_WP = find(droneData.data(:,2) == 999);
                idx_WP = [idx_start;idx_WP];
                c = 0:length(idx_WP)-1;
                idx_WP_corrected = idx_WP - c';   %to remove the shift caused by the 999
                
                % Find samples to ignore (sensor jumps), each pair indicates start and stop of
                % segment to ignore; remove all WP segments containing these indices
                WP_to_ignore = [];
                for w = 1:length(idxToRemove)/2
                    currIdxToRemove = idxToRemove((2*w-1):(2*w));
                    
                    WP_to_ignore_first = find(idx_WP_corrected<currIdxToRemove(1) ,1,'last');
                    if isempty(WP_to_ignore_first); WP_to_ignore_first = 1;end
                    
                    WP_to_ignore_last = find(idx_WP_corrected>currIdxToRemove(2) ,1,'first');
                    if isempty(WP_to_ignore_last); WP_to_ignore_last = 1;end
                    
                    WP_to_ignore = [WP_to_ignore,WP_to_ignore_first:WP_to_ignore_last];
                    
                end
                
                
                % initialize vectors
                [headRotAmp, torsoRotAmp, torsoSpeed_mean,torsoSpeed_max] = deal(nan(nWP-1,3));
                [headTorsoCorr,cc_peak] = deal(nan(nWP-1,5));
                [torsoSAL, headSAL, droneSAL, AI_roll, AI_pitch, AI_yaw,...
                    headPks_p,headPks_r, headPks_y, torsoPks_p, torsoPks_r, torsoPks_y,...
                    dronePks_r,dronePks_p,dronePks_y,torsoCorr_rp,torsoCorr_py,torsoCorr_ry,...
                    torsoSpeed_norm_mean, torsoSpeed_norm_max,...
                    dtw_roll, dtw_pitch, dtw_yaw] = deal(nan(nWP-1,1));
                
                nWPtoTake = min(length(idx_WP), nWP);
                % for drone metrics, use all WP
                for wp = 1:nWPtoTake-1
                    seg_idx = idx_WP_corrected(wp):idx_WP_corrected(wp+1);
                    
                    if length(seg_idx) > 1
                        %SAL
                        RPYrates = droneVals.data(seg_idx, 5:7);
                        droneSAL(wp,:) = sal3d( RPYrates(:,1)', RPYrates(:,2)', RPYrates(:,3)', SR);
                        
                        %Velocity peaks
                        dronePks_r(wp,:) = length(findpeaks(RPYrates(:,1)))+length(findpeaks(-RPYrates(:,1)));
                        dronePks_p(wp,:) = length(findpeaks(RPYrates(:,2)))+length(findpeaks(-RPYrates(:,2)));
                        dronePks_y(wp,:) = length(findpeaks(RPYrates(:,3)))+length(findpeaks(-RPYrates(:,3)));
                    end
                end
                
                % normalize number of peaks to time
                dronePks_r = dronePks_r./time(2:end);
                dronePks_p = dronePks_p./time(2:end);
                dronePks_y = dronePks_y./time(2:end);
                
                % s029 got lost during the trial -> non-representative
                % number of peaks, exclude
                if strcmp(subjects{ii},'s029')
                    dronePks_r(15) = nan;
                    dronePks_p(15) = nan;
                    dronePks_y(15) = nan;
                end
                
                % for body metrics, use only good WP
                good_wp = setdiff(1:nWPtoTake-1,WP_to_ignore);
                for ll = 1:length(good_wp)-1
                    wp = good_wp(ll);
                    
                    if good_wp(ll+1)-good_wp(ll) == 1 % only for consecutive wp
                        
                        seg_idx = idx_WP_corrected(wp):idx_WP_corrected(wp+1);
                        
                        if length(seg_idx) > 1
                            
                            %%%%%%% Correlation %%%%%%%%%%
                            % Correspondance: head2 = body3; head3 = body 1; head1 = body2; = [6;7;2]
                            % = [pitch; yaw; roll ]
                            % Order: h1b1 h1b2 h1b3 h2b1 h2b2 h2b3 h3b1 h3b2 h3b3
                            headRot_segment = bodyAngles.data(seg_idx, 1:3);     % X - Y - Z = [pitch; yaw; roll ]
                            bodyRot_segment= bodyAngles.data(seg_idx, 4:6);    % Roll - Pitch - Yaw
                            htCorr= reshape(corr(headRot_segment, bodyRot_segment)', 1,9);
                            headTorsoCorr(wp,:) = htCorr([2,4,6,7,9]) ;    % 1 = H pitch - B pitch, 2 = H roll - Body roll, 3 = H roll - B yaw, 4 = H yaw - B roll, 5 = H yaw - B yaw
                            
                            % Movement amplitude: quartiles
                            headRotAmp(wp,:) = iqr(headRot_segment); %1 = X = Pitch, 2 = Y = Roll, 3 = Z = yaw
                            torsoRotAmp(wp,:) = iqr(bodyRot_segment); % Roll - Pitch - Yaw
                            
                            % SAL
                            torsoSAL(wp,:) = sal3d(cent_diff_3(bodyRot_segment(:,1),SR)', cent_diff_3(bodyRot_segment(:,2),SR)',cent_diff_3(bodyRot_segment(:,3),SR)',SR);
                            headSAL(wp,:) = sal3d(cent_diff_3(headRot_segment(:,1),SR)', cent_diff_3(headRot_segment(:,2),SR)',cent_diff_3(headRot_segment(:,3),SR)',SR);
                            
                            
                            % Anchoring index
                            AI_roll(wp,:) = anchoringIndex(headRot_segment(:,2), bodyRot_segment(:,1));
                            AI_pitch(wp,:) = anchoringIndex(headRot_segment(:,1), bodyRot_segment(:,2));
                            AI_yaw(wp,:) = anchoringIndex(headRot_segment(:,3), bodyRot_segment(:,3));
                            
                            % Cross-correlation peak time
                            % if < 0 -> head leading
                            cc_idx_head =[1 2 3 2 3]; % p-p, r-r, y-y, r-y, y-r
                            cc_idx_torso = [2 1 3 3 1];
                            
                            l = min(15, floor(length(seg_idx)/2));
                            for c = 1:5
                                if c ==1
                                    [r, lags] = xcorr(headRot_segment(:,cc_idx_head(c)), -bodyRot_segment(:,cc_idx_torso(c)),l); % signs inverted for pitch
                                else
                                    [r, lags] = xcorr(headRot_segment(:,cc_idx_head(c)), bodyRot_segment(:,cc_idx_torso(c)),l);
                                end
                                cc_peak(wp,c) = lags(find(r == max(r)))/SR;
                            end
                            
                            % velocity peaks
                            headSpeed = cent_diff_3(headRot_segment,1/SR);
                            torsoSpeed = cent_diff_3(bodyRot_segment,1/SR);
                            headPks_p(wp,:) = length(findpeaks(headSpeed(:,1)))+length(findpeaks(-headSpeed(:,1)));
                            headPks_r(wp,:) = length(findpeaks(headSpeed(:,2)))+length(findpeaks(-headSpeed(:,1)));
                            headPks_y(wp,:) = length(findpeaks(headSpeed(:,3)))+length(findpeaks(-headSpeed(:,3)));
                            torsoPks_r(wp,:) = length(findpeaks(torsoSpeed(:,1)))+length(findpeaks(-torsoSpeed(:,1)));
                            torsoPks_p(wp,:) = length(findpeaks(torsoSpeed(:,2)))+length(findpeaks(-torsoSpeed(:,2)));
                            torsoPks_y(wp,:) = length(findpeaks(torsoSpeed(:,3)))+length(findpeaks(-torsoSpeed(:,3)));
                            
                            % correlation Torso rotations
                            torsoCorr_rp(wp) = corr(bodyRot_segment(:,1), bodyRot_segment(:,2));
                            torsoCorr_ry(wp) = corr(bodyRot_segment(:,1), bodyRot_segment(:,3));
                            torsoCorr_py(wp) = corr(bodyRot_segment(:,2), bodyRot_segment(:,3));
                            
                            % speeds
                            torsoSpeed_mean(wp,:) = mean(abs(torsoSpeed),'omitnan');
                            torsoSpeed_max(wp,:) = max(abs(torsoSpeed));
                            torsoSpeed_norm_mean(wp) = mean(sqrt(torsoSpeed(:,1).^2+torsoSpeed(:,2).^2+torsoSpeed(:,3).^2));
                            torsoSpeed_norm_max(wp) = max(sqrt(torsoSpeed(:,1).^2+torsoSpeed(:,2).^2+torsoSpeed(:,3).^2));
                            
                            %dtw: interpolate so that each segment has 100 points; zscore
                            % data
                            L = 100;
                            
                            t = 1:length(headRot_segment);
                            t_interp = linspace(1,length(headRot_segment),L);
                            
                            headVals_interp = interp1(t,headRot_segment,t_interp);
                            bodyVals_interp = interp1(t,bodyRot_segment,t_interp);
                            
                            dtw_roll(wp,:) = dtw(zscore(bodyVals_interp(:,1)),zscore(headVals_interp(:,2)));
                            dtw_pitch(wp,:) = dtw(zscore(bodyVals_interp(:,2)),zscore(headVals_interp(:,1)));
                            dtw_yaw(wp,:) = dtw(zscore(bodyVals_interp(:,3)),zscore(headVals_interp(:,3)));
                            
                        end
                    end
                end
                
                %normalize peaks to time
                headPks_p = headPks_p./time(2:end);
                headPks_r = headPks_r./time(2:end);
                headPks_y = headPks_y./time(2:end);
                torsoPks_p = torsoPks_p./time(2:end);
                torsoPks_r = torsoPks_r./time(2:end);
                torsoPks_y = torsoPks_y./time(2:end);
                
                % to avoid wrong (large) n peaks, ignore segments with time < 1
                idx_t00 = find(time(2:end) < 1);
                [torsoPks_y(idx_t00), torsoPks_r(idx_t00),torsoPks_p(idx_t00),...
                    headPks_y(idx_t00), headPks_r(idx_t00),headPks_p(idx_t00),...
                    dronePks_y(idx_t00), dronePks_r(idx_t00),dronePks_p(idx_t00)] = deal(nan);
                
                % Compute path efficiency
                % first interpolate to create ideal path with Catmull-Rom splines
                
                % Add starting position and additional points before and after
                % trajectory, for interpolation
                WP0 = droneVals.data(1,8:10); % drone starting position
                WP00 = WP0 - (WPInfo.data(1,:) - WP0);
                WP_end = WPInfo.data(end,:) + ( WPInfo.data(end,:) - WPInfo.data(end-1,:));
                WP = [WP00; WP0; WPInfo.data;WP_end];
                
                Tension = 0.5; % curvature ratio
                n = 75; % n points
                idealPath = [];
                for wp = 2:length(WP)-2
                    P0 = WP(wp-1,:);
                    P1 = WP(wp,:);
                    P2 = WP(wp+1,:);
                    P3 = WP(wp+2,:);
                    
                    idealPath_segment=crdatnplusoneval(P0,P1,P2,P3,Tension,n)';
                    idealPath = [idealPath;idealPath_segment];
                end
                
                dx = diff(droneVals.data(:,8:10));
                path_length_real = sum(sqrt(sum(dx.^2,2)));
                dx_ideal = diff(idealPath);
                path_length_ideal = sum(sqrt(sum(dx_ideal.^2,2)));
                
                path_ratio = path_length_real/path_length_ideal;
                
                
                % Add to table...
                T_ ={subjects{ii},num2str(ages(ii)),genders(ii), (cellstr(control_mode(kk))),(cellstr(phases(p))),dist2WP,time,...
                    -headTorsoCorr(:,1),headTorsoCorr(:,2),headTorsoCorr(:,3),headTorsoCorr(:,4),headTorsoCorr(:,5),...
                    headRotAmp(:,1),headRotAmp(:,2),headRotAmp(:,3),torsoRotAmp(:,1),torsoRotAmp(:,2),torsoRotAmp(:,3),...
                    torsoSAL, headSAL, droneSAL, AI_roll, AI_pitch, AI_yaw,...
                    cc_peak(:,1), cc_peak(:,2), cc_peak(:,3), cc_peak(:,4), cc_peak(:,5),...
                    headPks_p, headPks_r, headPks_y, torsoPks_r,torsoPks_p, torsoPks_y,...
                    dronePks_r,dronePks_p,dronePks_y,torsoCorr_rp,torsoCorr_py,torsoCorr_ry,...
                    torsoSpeed_mean(:,1), torsoSpeed_mean(:,2), torsoSpeed_mean(:,3), ...
                    torsoSpeed_max(:,1), torsoSpeed_max(:,2), torsoSpeed_max(:,3),...
                    torsoSpeed_norm_mean, torsoSpeed_norm_max,errorDir(:,1), errorDir(:,2),...
                    dtw_roll, dtw_pitch, dtw_yaw,path_ratio};
                
                T = [T;T_];
                
                clear head* torso* *SAL AI*
            end
        end
    end
    
    
end


save([saveFolder,'/dataTable.mat'],'T');


%% Average values for each field by session
T_mean = T(:,1:5);
VarNames = T.Properties.VariableNames;
for ii = 6:width(T)-1   % last field is path_efficiency
    T_mean.(VarNames{ii}) = cellfun(@(x)mean(x,'omitnan'),T.(VarNames{ii}));
end
T_mean.Path_ratio = T.Path_ratio;

save([saveFolder,'/dataTable_means.mat'],'T_mean');



