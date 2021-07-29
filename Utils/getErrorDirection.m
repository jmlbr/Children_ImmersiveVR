
function errorDir = getErrorDirection(WPInfo, droneVals, dist2WP, nWP)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the 2D error direction for each waypoint
% Plot :   compass(errorDir(:,1), errorDir(:,2))
% On the plot, the arrows indicate the direction in which the drone SHOULD
% go (pointing downwards = too high; pointing to the left = too much on the
% right)
% Arguments:
% WPInfo: importdata(xxx_manoeuvreList.txt)
% droneVals: readDroneAngles(droneFile.txt)
% dist2WP: distance to each Waypoint
% nWP: number of waypoints; added separately since sometimes Unity stacks 2
% series together...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idxWP = 1:min(nWP,length(WPInfo.data))
    dist2currWP = sqrt(sum((repmat(WPInfo.data(idxWP,:),length(droneVals.data),1)-droneVals.data(:,8:10)).^2,2));
    
    [~,idxMinDist] = min(abs(dist2currWP - dist2WP(idxWP)));   % find when distance to WP is closest to value computed by unity when plane is crossed
    
    direction = WPInfo.data(idxWP,:)-droneVals.data(idxMinDist,8:10);
    
    % find sidewards direction of error
    % sign > 0 means drone is on the right of the path (should go to the
    % left)
    if idxWP > 1
      v1 = [WPInfo.data(idxWP,[1,3])-WPInfo.data(idxWP-1, [1,3]),0];
      v2 = [droneVals.data(idxMinDist,[8,10])-WPInfo.data(idxWP-1, [1,3]),0];
      
      crossprod = cross(v2, v1);
      signs(idxWP) = sign(crossprod(3));
      lateral_error = sign(crossprod(3))*sqrt(direction(1)^2+direction(3)^2);
      
    else
      lateral_error = 0;
    end
    
    % vertical error > 0 means drone is below path (should go up)
    errorDir(idxWP,:) = [lateral_error, direction(2)];
end
  
% 
% figure;
%   compass(errorDir(:,1), errorDir(:,2))
%   hold on
%   h = compass(median(errorDir(:,1)), median(errorDir(:,2)),'r');
%   h.LineWidth = 2;
%   
  
  