clear
%% Read h5 files. This function is in /Users/jun/Documents/Work/Project/Licking/Analysis_Code
[areas, Jaw_heights] = extract_h5();
%% Plot traces for inspection
% use (smoothdata(areas,'gaussian', 12)) for plotting small time window. if
% no parameter is set, this command makes trace too smooth. Or consider to
% use different kind of filters.

plot(smoothdata(areas,'gaussian', 12))
hold on
plot(5*smoothdata(Jaw_heights{1},'gaussian', 12))
ylim([0 15000])
%% Store data in a structure file

S.areas = areas;
S.Jaw_heights = Jaw_heights{1};
%% Concat data - skip if only one set of h5 files is analyzed
% Repeat this section to add h5 files

% Perhaps the most accurate way to store data for different videos is to 
% assign a specific row to each video. However, in this case, concatenating
% the data vertically would not affect the accuracy.

clearvars -except S
[areas, Jaw_heights] = extract_h5();

S.areas = vertcat(S.areas, areas);
S.Jaw_heights = vertcat(S.Jaw_heights, Jaw_heights{1});
%% use this if start from an S file
areas = S.areas;
Jaw_heights = S.Jaw_heights;
%%
[peaks,locs] = findpeaks(smoothdata(areas,'gaussian', 12),'MinPeakDistance',17.5, 'MinPeakHeight',200);

% MinPeakDistance @350fps. 2.85714 ms/frame. Cut off licking happens faster
% than 20Hz. Over @20Hz, the n + 1 lick happens within 50ms after n lick. 50/2.85714 = 17.5
% Licking happens this MinPeakDistance must be artifact.

figure;
findpeaks(smoothdata(areas,'gaussian', 12),'MinPeakDistance',17.5,'MinPeakHeight',200);
% xlim([0,frames(end)]);
xlim([0,length(areas)]);
ylim([0,15000]);
xlabel('Frame');
ylabel('Area');

%%
% time window for consecutive spikes. 70 frames = 200ms corresponds to 5Hz 

% Peak locations (in frames)
peak_locs = locs;

% maximum distance between peaks to be considered in the same aggregate
max_distance = 70;

% initialize variables
agg_count = 0;
peak_count = 0;
agg_locs = [];
num_peaks = [];

% loop through each peak location
for i = 1:length(peak_locs)
    
    % if this is the first peak, initialize the aggregate
    if peak_count == 0
        agg_count = agg_count + 1;
        agg_locs{agg_count} = peak_locs(i);
        peak_count = peak_count + 1;
        
    % if this peak is within the maximum distance of the previous peak,
    % add it to the same aggregate
    elseif peak_locs(i) - agg_locs{agg_count}(end) <= max_distance
        agg_locs{agg_count} = [agg_locs{agg_count} peak_locs(i)];
        peak_count = peak_count + 1;
        
    % if this peak is too far from the previous peak, start a new aggregate
    else
        % store the number of peaks in the previous aggregate
        num_peaks(agg_count) = peak_count;
        
        agg_count = agg_count + 1;
        agg_locs{agg_count} = peak_locs(i);
        peak_count = 1;
    end
    
end

% store the number of peaks in the final aggregate
num_peaks(agg_count) = peak_count;

% print the aggregate locations and number of peaks
for i = 1:length(agg_locs)
    fprintf('Aggregate %d contains %d peaks: %s\n', ...
        i, num_peaks(i), mat2str(agg_locs{i}));
end

LickPerBout = mean(num_peaks);
%% Inter Lick Interval
InterLick_interval = diff(peak_locs);
InterLick_interval = InterLick_interval(InterLick_interval <= 700); % cut off 2 second
InterLick_interval_InSecond = InterLick_interval/350;
InterLick_interval_InMillisecond = InterLick_interval_InSecond*1000;

[histFreq, histXout] = hist(InterLick_interval_InMillisecond, 100); % 20 ms bin 
figure;
bar(histXout, histFreq/sum(histFreq)*100);
xlabel('Inter Lick Interval (ms)');
ylabel('Probability');

xlim([0,2000]);
ylim([0,50]);
%% Find the right time window for tresholding
% close()
% [yupper,ylower] = envelope(smoothdata(areas,'gaussian'),175,'rms');
% plot(smoothdata(areas,'gaussian'))
% hold on
% plot(yupper)
%% Detect licking period

% The size of the detection windows varies with the size of the time window (set at 175 for now). 
% The start and end of the window may need to be adjusted Using the information from the first and last licks. 
% For example, the window begins 200ms before the peak amplitude of the first lick.

[yupper,ylower] = envelope(smoothdata(areas,'gaussian'),175,'rms'); 
lowthreshold = min(yupper) + min(yupper)*0.2; %This part has to be adaptive. May need to change.
aboveThreshold = (yupper > lowthreshold);
aboveThreshold = [false, aboveThreshold', false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);  
spanWidth = falling - rising;
wideEnough = spanWidth >= 70; % longer than 0.2 sec   
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
%% Store licking and chewing period in cells
licking_period = cell(length(startPos),1);
chewing_period = cell(length(startPos),1);

for i = 1:length(startPos)
    licking_period{i} = areas(startPos(i):endPos(i));
    chewing_period{i} = Jaw_heights{1}(startPos(i):endPos(i));
end

%% Pick licking period that contain more than 3 licks 

num_peaks = zeros(length(startPos),1);
for i = 1:length(startPos)
    [peaks_all,locs_all] = findpeaks(smoothdata(licking_period{i},'gaussian'),'MinPeakDistance',20, 'MinPeakHeight',lowthreshold);
    num_peaks(i) = numel(peaks_all);
    LickPeaks{i} = peaks_all; 
end

ind = find(num_peaks >= 3);
LickPeaksConcat = cat(1, LickPeaks{:});

licking_multi = {licking_period{ind,:}};
chewing_multi =  {chewing_period{ind,:}};
%% For visual inspection of the detected licking periods
k = 18;
plot(smoothdata(licking_period{k},'gaussian'))
hold on
plot(5*smoothdata(chewing_period{k},'gaussian'))
findpeaks(smoothdata(licking_period{k},'gaussian'),'MinPeakDistance',25, 'MinPeakHeight',lowthreshold);
%% Correlation coefficient between tongue and jaw
corr_lick_chew = zeros(length(licking_multi),1);

for i = 1: length(licking_multi)
    R = corrcoef(smoothdata(licking_multi{i},'gaussian'), smoothdata(chewing_multi{i},'gaussian'));
    corr_lick_chew(i) = R(1,2);
end
Mean_corr = mean(corr_lick_chew);

%% Licking frequency
for i = 1 : length(licking_multi)
    threshold = min(licking_period{i});
    aboveThreshold = (smoothdata(licking_multi{i},'gaussian') > threshold);
    aboveThreshold = [false, aboveThreshold', false];  %pad with 0's at ends
    edges = diff(aboveThreshold);
    first = find(edges, 1, 'first');
    last= find(edges, 1, 'last');
    duration = (last - first)/350;
    [peaks,locs] = findpeaks(smoothdata(licking_multi{i},'gaussian'),'MinPeakDistance',25, 'MinPeakHeight',lowthreshold);
    NumberOfLicks = length(peaks);
    LickPerBout(i) = length(peaks); 
    LickFrequency(i) = NumberOfLicks/duration;
end
Mean_LickFrequency = mean(LickFrequency);
Mean_LickPerBout = mean(LickPerBout);
%% jaw amplitude
chewing = cat(1, chewing_multi{:}); % Setpoint detection deosn't work well wwithh small pieces of chewing period. Concat all the chewing period. 

Fs = 350;
setpt_cut = 1; % minimum frequency
lowpass_cut = 15; % maximum frequency 40 Hz (for rat and mouse)
filter_order = 3;

[bh, ah] = butter(filter_order,setpt_cut/(Fs/2),'low');
setpoint = filtfilt(bh,ah, smoothdata(chewing,'gaussian'));
chewing_setpoint = smoothdata(chewing) - setpoint;
plot(chewing)
hold on
plot(chewing_setpoint)
plot(setpoint)

amp_pp = fan_whiskutil_get_whiskamplitude_DH(chewing_setpoint);
avg_amp = mean(amp_pp);
%% Store data

S.Peaks = LickPeaks;
S.PeaksConcat = LickPeaksConcat;
S.Corr_lick_chew = corr_lick_chew;
S.Mean_Corr_lick_chew = Mean_corr;
S.LickFrequency = LickFrequency;
S.Mean_LickFeaquency = Mean_LickFrequency; 
S.LickPerBout = LickPerBout;
S.Mean_LickPerBout = Mean_LickPerBout;
S.Jaw_amplitude = avg_amp;

filename = fullfile(file_tongue);
filename = filename(1:end-6);
filename = strcat(filename, '_kinematics');
save(filename,'S')
%%



%% load kinematics 1. Open the first structue file and run this
Peaks = S.PeaksConcat;
Corr_lick_chew = S.Corr_lick_chew;
Mean_Corr_lick_chew = S.Mean_Corr_lick_chew;
LickFrequency = S.LickFrequency;
Mean_LickFrequency = S.Mean_LickFeaquency;
LickPerBout = S.LickPerBout;
Mean_LickPerBout = S.Mean_LickPerBout;
Jaw_amplitude = S.Jaw_amplitude;
clear S
%% load kinematics 2. Open the second structure file and rund this. Then go the third run and so on.
Peaks = vertcat(Peaks, S.PeaksConcat);
Corr_lick_chew = vertcat(Corr_lick_chew, S.Corr_lick_chew);
Mean_Corr_lick_chew = vertcat(Mean_Corr_lick_chew, S.Corr_lick_chew);
LickFrequency = horzcat(LickFrequency, S.LickFrequency);
Mean_LickFrequency = vertcat(Mean_LickFrequency, S.Mean_LickFeaquency);
LickPerBout = horzcat(LickPerBout, S.LickPerBout);
Mean_LickPerBout = vertcat(Mean_LickPerBout, S.Mean_LickPerBout);
Jaw_amplitude = vertcat(Jaw_amplitude, S.Jaw_amplitude);
clear S
%% 
%% Before
Before.note = "Merge: 7/20,7/21";
Before.Peaks = Peaks;
Before.Corr_lick_chew = Corr_lick_chew;
Before.Mean_Corr_lick_chew = mean(Mean_Corr_lick_chew);
Before.LickFrequency = LickFrequency;
Before.Mean_LickFrequency = mean(Mean_LickFrequency);
Before.LickPerBout = LickPerBout;
Before.Mean_LickPerBout = mean(Mean_LickPerBout);
Before.Jaw_amplitude = mean(Jaw_amplitude);

%% One week
OneWeek.note = "Merge: 7/29";
OneWeek.Peaks = Peaks;
OneWeek.Corr_lick_chew = Corr_lick_chew;
OneWeek.Mean_Corr_lick_chew = mean(Mean_Corr_lick_chew);
OneWeek.LickFrequency = LickFrequency;
OneWeek.Mean_LickFrequency = mean(Mean_LickFrequency);
OneWeek.LickPerBout = LickPerBout;
OneWeek.Mean_LickPerBout = mean(Mean_LickPerBout);
OneWeek.Jaw_amplitude = mean(Jaw_amplitude);

%% After
After.note = "Merge: 8/8, 8/10";
After.Peaks = Peaks;
After.Corr_lick_chew = Corr_lick_chew;
After.Mean_Corr_lick_chew = mean(Mean_Corr_lick_chew);
After.LickFrequency = LickFrequency;
After.Mean_LickFrequency = mean(Mean_LickFrequency);
After.LickPerBout = LickPerBout;
After.Mean_LickPerBout = mean(Mean_LickPerBout);
After.Jaw_amplitude = mean(Jaw_amplitude);
