% This script is for analyzing phase relationship between chewing ans
% licking

TongueAreas = S.areas;
JawHeights = S.Jaw_heights;

JawHeights = double(JawHeights);
%% Apply Filters to Jaw Trace
samplingRate = 350;
setpointCutoff = 1;  % Minimum frequency
lowpassCutoff = 15;  % Maximum frequency 
filterOrder = 3;

% Setpoint
[bh, ah] = butter(filterOrder, setpointCutoff / (samplingRate / 2), 'low');
JawHeightsSetpoint = filtfilt(bh, ah, JawHeights);

% Lowpass 
[bl, al] = butter(filterOrder, lowpassCutoff / (samplingRate / 2), 'low');
JawHeightsLowpass = filtfilt(bl, al, JawHeights);

JawHeightsFiltered = (JawHeightsLowpass - JawHeightsSetpoint)';

%% Apply Filters to Tongue Trace
% Removing setpoint from licking trace seems to introduce artifact to the trace.
% Applying a lowpass filter creates artificial dips during the unseen licking period. 
% These dips do not appear to introduce strange artifact.

% Lowpass 
[bl, al] = butter(filterOrder, lowpassCutoff / (samplingRate / 2), 'low');
TongueAreasFiltered = filtfilt(bl, al, TongueAreas);

%% Plot to Test Filtering
plot(JawHeightsFiltered)
hold on
plot(TongueAreasFiltered*0.02)

%% Phase
% This part works fine. Just need to apply this during the licking period.

ChewingPhase = angle(hilbert(JawHeightsFiltered));
LickingPhase = angle(hilbert(TongueAreasFiltered));

plot(ChewingPhase)
hold on
plot(LickingPhase)