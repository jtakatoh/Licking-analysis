clear
%% Read h5 files.
image_height = 480;
image_width = 640;

% the file ends with 0 is tongue 
[file_tongue,path_tongue] = uigetfile('*0.h5');
Tongue = fullfile(path_tongue, file_tongue);
% Tongue_filename = '/Users/jun/Documents/Work/Project/Licking/Analysis_Code/Tongue-Jaw 2022_0902/Phox2B#8_20220720_1_1_0.h5';
% the file ends with 1 is jaw 
[file_jaw,path_jaw] = uigetfile('*1.h5');
Jaw = fullfile(path_jaw, file_jaw);
% Jaw_filename = '/Users/jun/Documents/Work/Project/Licking/Analysis_Code/Tongue-Jaw 2022_0902/Phox2B#8_20220720_1_1_1.h5';

frames=h5read(Tongue,'/frames');
heights=h5read(Tongue,'/heights');
widths=h5read(Tongue,'/widths');
probs=h5read(Tongue,'/probs');

Jaw_heights=h5read(Jaw,'/heights');


start_frame = 1;
end_frame = frames(end);
total_frames = end_frame - start_frame;

areas = zeros(end_frame - start_frame,1);
averages = zeros(end_frame - start_frame,2);

img = zeros(image_height,image_width);
img_array = cell(length(frames),1);

for j=1:length(frames)
    img(:) = 0;
    frame_ind = j; %This is the index of the frame output from the neural network
    %The actual frame number of the entire video is frames{frame_ind}

    for i=1:length(heights{frame_ind})
        if heights{frame_ind}(i) == 0
            heights{frame_ind}(i) = 1;
        end
        if widths{frame_ind}(i) == 0
            widths{frame_ind}(i) = 1;
        end
        img(heights{frame_ind}(i),widths{frame_ind}(i)) = probs{frame_ind}(i);
    end
    
    img(img<0.2) = 0; % remove low probability areas

    
    img_array(j) = {img};
    areas(frames(j),1) = sum(sum(img));
end

%%
plot(smoothdata(areas(5000:8500),'gaussian',7))
hold on
plot(5*smoothdata(Jaw_heights{1}(5000:8500),'gaussian',7))
ylim([0 15000])
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



%% load kinematics 1 
Peaks = S.PeaksConcat;
Corr_lick_chew = S.Corr_lick_chew;
Mean_Corr_lick_chew = S.Mean_Corr_lick_chew;
LickFrequency = S.LickFrequency;
Mean_LickFrequency = S.Mean_LickFeaquency;
LickPerBout = S.LickPerBout;
Mean_LickPerBout = S.Mean_LickPerBout;
Jaw_amplitude = S.Jaw_amplitude;
clear S
%% load kinematics 2
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



