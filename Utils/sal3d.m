function S = sal3d( vx, vy, vz, Ts, parameters )
% SPECTRALARCLENGTH computes the smoothness of the given movement speed 
% profile using the spectral arc length method.
% This function takes three inputs and provides one output.
% Inputs: { speed, Ts*, parameters* } ('*' optional parameters)
%         vx, vy, vz : The speed profile of the movement in x-,y, and z-axis. This is a 1xN
%         row vector. N is the total number of points in the speed profile.
%         The function assumes that the movement speed profile is already
%         filtered and segmented.
%
%         TS*: Sampling time in seconds. (DEFAULT VALUE = 100Hz) NOTE: IF
%         YOUR DATA WAS NOT SAMPLED AT 100HZ, YOU MUST ENTER THE
%         APPROPRIATE SAMPLING TIME FOR ACCURATE RESULTS.
%
%         PARAMETERS*: This contains the parameters to be used spectral arc
%         length computation. This is a 1x2 column vector. This input
%         argument is option
%           - PARAMETER(1): Cut-off frequency for the spectral arc length
%           calculation. (DEFAULT VALUE = 20HZ) NOTE: 20Hz IS USED TO
%           REPRESENT THE MAXIMUM FREQUENCY COMPONENT OF A MOVEMENT. THIS
%           WAS CHOSEN TO COVER BOTH NORMAL AND ABNORMAL MOTOR BEHAVIOUR.
%           YOU CAN USE A VALUE LOWER THAN 20Hz IF YOU ARE AWARE OF THE
%           MAXIMUM FREQUENCY COMPONENT IN THE MOVEMENT OF INTEREST IS
%           LOWER THAN 20 Hz.
%           - PARAMETER(2): Zero padding index. This parameter controls the
%           resolution of the movement spectrum calculated from the speed
%           profile. (DEFAULT VALUE = 4). NOTE: IT IS NOT ADVISABLE TO USE
%           VALUES LOWER THAN 4.
%
% Outputs: { S }
%          S: This is smoothness of the given movement.
%
% For any queries about the method or the code, or if you come across any 
% bugs in the code, feel free to contact me at siva82kb@gmail.com
% Sivakumar Balasubramanian. July 22, 2011. 

% Calculate speed vector
speed = sqrt(vx.^2+vy.^2+vz.^2);

% Check input arguments.
if nargin < 3
    disp('Error! Input at least the movement speed profile for the smoothness calculation.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
elseif nargin == 3
    % Default sampling time.
    Ts = 1/1000; % 1ms.
elseif nargin == 4
    % Default parameters are use for the spectral arc length caclulations.
    parameters = [20,4];    
end;

% Check if the input argument are of the appropriate dimensions.
% Speed profile.
if ~isrow(speed)
    disp('Error! speed must be a row vector.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;
% Sampling time.
if ~isscalar(Ts)
    disp('Error! Ts must be a scalar.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;
% Parameters.
if length(parameters) ~= 2 
    disp('Error! parameter is a vector with two elements.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;

% Calculate the spectrum of the speed profile.
N = length(speed);
Nfft = 2^(ceil(log2(N))+parameters(2));
speedSpectrum = abs(fft( speed, Nfft ));

% Normalize spectrum with respect to the DC component.
speedSpectrum = speedSpectrum/speedSpectrum(1);

% Get index corresponding to the cut off frequency.
freq = 0:(1/Ts)*(1/Nfft):(1/Ts)*((Nfft-1)/Nfft);
inxFc = find( freq(1:end) <= parameters(1), 1, 'last' );

% Calculate the spectral arc length.
% 1. select the spectrum of interest.
speedSpectrum = speedSpectrum(1:inxFc);
% 2. Calculate the incremental arc lengths.
dArcLengths = sqrt((1/(inxFc-1))^2 + (diff(speedSpectrum)).^2);
% 4. Compute movement smoothness.
S = -sum(dArcLengths);

return;