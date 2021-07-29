% createDataTable_JAR.m

% Processing routine to extract kinematic parameters from the files
% recorded during the JAR experiment.
% Variables:
% - Average Error
% - Variability of the error
% - Head-torso correlation and DTW distance (roll plane)
% - Head anchoring index (roll plane)
% - Head-torso amplitude difference (roll plane)
% - Average overshoot
% - Number of oscillations around target value
% - "Error" of the non-controlling body part


% Saves one data table containing average values for each block

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Head-torso plane correspondance:
% torso Roll = torso x, head roll = head y

% ©Jenifer Miehlbradt, EPFL, 2021

%% Modify paths here
clearvars
clc

% Helper functions
addpath('Utils')   % ADD 'Utils' FOLDER HERE;


mainFolder = cd 
dataFolder = [mainFolder,'/rawFiles/']; 
saveFolder = cd; % ADD FOLDER TO SAVE DATA TABLES HERE;
cd(dataFolder)

%% Get unique subject identifiers
bodyAngleFiles = dir([dataFolder,'S*BodyAngles*']);
bodyAngleFiles = {bodyAngleFiles.name};

fileNameParts = cellfun(@(x) strsplit(x, '_'), bodyAngleFiles, 'UniformOutput', false);
fileNameParts = vertcat(fileNameParts{:});
[~,idx]= (unique(cell2mat(fileNameParts(:,1)),'rows'));
subjects = fileNameParts(idx,1);
clear idx fileNameParts


%% Initialize
% Load subject information
load([mainFolder,'/subjectInfo.mat']);

% Conditions
blocks = {'training_head_still','training_torso_still', 'test_head_still','test_torso_still', 'test_head_forward','test_torso_forward'};

% Window in which the final, stable position is calculated
window = 22; %last 1.5 sec of the trial

% Sampling rate
SR = 1000/68;    % 1 sample all 68 ms =14.7 Hz

% 
[ID, age, gender, control, session, ...
    avgError, varError,...
    HTCorrRoll, DTWRoll, AIRoll,...
    ampDiffRoll,avgOvershoot,...
    avgOscillations, ...
    ncError] = deal([]);

T_JAR = table(ID, age, gender, control, session, ...
    avgError,varError,...
    HTCorrRoll, DTWRoll, AIRoll,...
    ampDiffRoll,avgOvershoot,...
    avgOscillations,...
    ncError,...
    'VariableNames',{'ID', 'Age', 'Gender', 'Control', 'Session', ...
    'avgError', 'varError',...
    'HTCorrRoll', 'DTWRoll', 'AIRoll',...
    'ampDiffRoll','avgOvershoot',...
    'avgOscillations',...
    'ncError'});

%%
for s = 1:length(subjects)    
    disp(['Processing subject ', num2str(s)])
    for b = 1:length(blocks) % blocks
        subID = subjects{s};
        infoIdx = find(ismember(subjectInfo.ID,subID));
        common_name = [subID,'_', blocks{b}];
        angleFile = dir([common_name,'_BodyAngles*.txt']);
        eventFile = dir([common_name,'_Events*.txt']);
        control = [];
        
        if ~isempty(angleFile)
            if contains(common_name, 'training')
                session = 'feedback';
            elseif contains(common_name,'forward')
                session = 'forward';
            else
                session = 'still';
            end
            
            % Load data
            events = txtToStruct(eventFile.name);   % Trigger values with target angles
            angles = txtToStruct(angleFile.name);   % Head and body angles (xyz) and time vector
            
            fields = fieldnames(angles);
            if contains(common_name, 'head')
                controlData = angles.(fields{3});
                nonControlData = angles.(fields{5});
                control = 'Head';
            elseif contains(common_name, 'torso')
                controlData = angles.(fields{5});
                nonControlData = angles.(fields{3});
                control = 'Torso';
            end
            
            % Construct continuous vector of target angles
            targetVals = nan(1,length(angles.Time));
            for ii = 1:length(events.EventTime)-1
                idx_event(ii) = find (angles.Time > events.EventTime(ii),1,'first');
            end
            idx_event(end+1) = length(angles.Time);
            for ii = 1:length(events.EventTime)-1
                targetVals(idx_event(ii):idx_event(ii+1)) = events.TargetAngle(ii);
            end
            
            % Find data corresponding to each trial, compute error and
            % "error" of the non-controlling body part
            toTake = -ones(length(events.EventTime)-1,1);
            [trialAngle,trialError,nonControllingError,ampDiff,overshoot,...
                damping,duration,noscillations] = deal(nan(length(events.EventTime)-1,1));
            
            for ii = 1:length(events.EventTime)-1
                i1 = find (angles.Time > events.EventTime(ii),1,'first');
                i2 = find (angles.Time < events.EventTime(ii+1),1,'last');
                duration(ii) = (angles.Time(i2)-angles.Time(i1))/1000;
                eventData = controlData(i1:i2);
                ncEventData = nonControlData(i1:i2);
                
                if duration(ii) > 2.5 && duration(ii) < 8 ... % For early trials (w/o fixed duration), take trials lasting between 2.5 and 8 seconds
                        && all(abs(eventData) < 35)               % Exclude abnormal range
                    if (ii  > 1 && events.TargetAngle(ii)~=events.TargetAngle(ii-1)) ||...
                            (ii  == 1 && events.TargetAngle(ii)~=0)  % Exclude trials where end and start angle are the same
                        toTake(ii) = 1;
                        
                        % Align all trials towards the same direction
                        if (ii  > 1 && events.TargetAngle(ii)<events.TargetAngle(ii-1)) ||...
                                (ii  == 1 && events.TargetAngle(1)== -15)
                            eventData = -eventData;
                            ncEventData = -ncEventData;
                        end
                        
                        stableIdx = (length(eventData)-window):length(eventData);
                        trialAngle(ii) = mean(eventData(stableIdx));
                        trialError(ii) = trialAngle(ii)-abs(events.TargetAngle(ii));
                        if abs(trialError(ii)) < 15  % If error > 15, the participant moved towards the wrong direction
                            
                            % Compute "error" (alignment) of the non-controlling body part
                            nonControllingError(ii) = mean(ncEventData(stableIdx))-abs(events.TargetAngle(ii));
                            
                            % Compute amplitude difference between both body parts
                            if strcmp(control, 'Head')
                                ampDiff(ii) = mean(eventData(stableIdx))-mean(ncEventData(stableIdx));
                            elseif strcmp(control,'Torso')
                                ampDiff(ii) = mean(ncEventData(stableIdx))-mean(eventData(stableIdx));
                            end
                            
                            % Compute overshoot
                            overshoot(ii) = max(0,max(eventData) - trialAngle(ii));
                            % Compute number of oscillations = Number of crossings of the finalvalue / 2
                            noscillations(ii) = ceil(0.5*numel(find((eventData(1:end-1)-trialAngle(ii)).*(eventData(2:end)-trialAngle(ii)) < 0)));
                            
                        else   % If wrong direction, exclude trials
                            toTake(ii) = 0;
                        end
                    else
                        toTake(ii) = -3;
                    end
                end
            end
            
            % Compute, average and store metrics
            % Error
            varError = std(trialError,'omitnan');
            avgError = mean(trialError, 'omitnan');
            
            % Non-controlling error
            ncError = mean(nonControllingError,'omitnan');
            
            % Amp diff
            ampDiffRoll = mean(ampDiff,'omitnan');
            
            % Overshoot
            avgOvershoot = mean(overshoot,'omitnan');
            
            % Oscillations
            avgOscillations = mean(noscillations,'omitnan');
            
            % Head-Torso correlations
            idx_start = find (angles.Time > events.EventTime(1),1,'first');
            HTCorrRoll = corr(angles.(fields{3})(idx_start:end)', angles.(fields{5})(idx_start:end)');
            
            
            % Roll DTW
            headRoll_interp = interp1(1:length(angles.Time(idx_start:end)),angles.(fields{3})(idx_start:end),linspace(1,length(angles.Time(idx_start:end)),2000));
            torsoRoll_interp = interp1(1:length(angles.Time(idx_start:end)),angles.(fields{5})(idx_start:end),linspace(1,length(angles.Time(idx_start:end)),2000));
            DTWRoll= dtw(zscore(headRoll_interp),zscore(torsoRoll_interp));
            
            % Head AI Roll
            AIRoll = anchoringIndex(angles.(fields{3})(idx_start:end)', angles.(fields{5})(idx_start:end)');
            
            if s == 11 && b ==6 % missing data
                [HTCorrRoll,DTWRoll,AIRoll,nc] = deal(nan);     % Head data not recorded
            end
            
            
            t_ = {subID, subjectInfo.Age(infoIdx),subjectInfo.Gender(infoIdx),control,session,...
                avgError, varError,HTCorrRoll, DTWRoll, AIRoll,...
                ampDiffRoll,avgOvershoot,avgOscillations,ncError};
            
            T_JAR = [T_JAR;t_];
        end
    end
end

save([saveFolder,'/T_JAR.mat'],'T_JAR')


