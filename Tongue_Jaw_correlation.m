clear
%% Read h5 files.
image_height = 480;
image_width = 640;

% the file ends with 0 is tongue 
Tongue_filename = '/Users/jun/Documents/Work/Project/Licking/Analysis_Code/Tongue-Jaw 2022_0902/Phox2B#8_20220720_1_1_0.h5';
% the file ends with 1 is jaw 
Jaw_filename = '/Users/jun/Documents/Work/Project/Licking/Analysis_Code/Tongue-Jaw 2022_0902/Phox2B#8_20220720_1_1_1.h5';

frames=h5read(Tongue_filename,'/frames');
heights=h5read(Tongue_filename,'/heights');
widths=h5read(Tongue_filename,'/widths');
probs=h5read(Tongue_filename,'/probs');

Jaw_heights=h5read(Jaw_filename,'/heights');


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
lowthreshold = min(yupper) + min(yupper)*0.1; %This part has to be adaptive. May need to change.
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
    [peaks,locs] = findpeaks(smoothdata(licking_period{i},'gaussian'),'MinPeakDistance',35, 'MinPeakHeight',lowthreshold);
    num_peaks(i) = numel(peaks);
end

ind = find(num_peaks >= 3);

licking_multi = {licking_period{ind,:}};
chewing_multi =  {chewing_period{ind,:}};
%% For visual inspection of the detected licking periods
k = 1;
plot(smoothdata(licking_period{k},'gaussian'))
hold on
plot(5*smoothdata(chewing_period{k},'gaussian'))
% findpeaks(smoothdata(licking_period{k},'gaussian'),'MinPeakDistance',35, 'MinPeakHeight',lowthreshold);
%% Correlation coefficient between tongue and jaw
corr_lick_chew = zeros(length(licking_multi),1);

for i = 1: length(licking_multi)
    R = corrcoef(smoothdata(licking_multi{i},'gaussian'), smoothdata(chewing_multi{i},'gaussian'));
    corr_lick_chew(i) = R(1,2);
end
Mean_corr = mean(corr_lick_chew);
