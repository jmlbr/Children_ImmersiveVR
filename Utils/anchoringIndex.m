% Computes the anchoring index between 2 segments (torso = reference, 
% head = "moving segment")
% Ref: Assaiante 2012, PLoS ONE

function AI = anchoringIndex(headAngle, torsoAngle)
headAngle_rel = headAngle-torsoAngle;

std_rel = std(headAngle_rel);
std_abs = std(headAngle);

AI = (std_rel^2 - std_abs^2)/(std_rel^2 + std_abs^2);
end

