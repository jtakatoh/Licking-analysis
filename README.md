# Licking-analysis

## install ffmpeg
https://cyber-tenchou.com/tips/20200904/  
https://eng-entrance.com/linux-vi-save  
vi .zshrc  
Write: export PATH=/Applications:$PATH  
To save and exit: :wq 

## Making masks

## FFMPEG  

### Check the number of total frames  
ffprobe -v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 /home/wanglab/Desktop/build-CameraViewer-Qt_Static-Debug/Phox2b_#14_teLC_Choco_20220908_1_350fps0.mp4  

### Extract frames from the video (original or tracked/labeled video).  
VLC can extract frames but it drops frames. ffmpeg does a better job.  
Change directory to were you want to output frame images and run this (bmp is better than png)  
ffmpeg -r 1 -i /home/wanglab/Desktop/build-CameraViewer-Qt_Static-Debug/Phox2b_#14_teLC_Choco_20220908_3_350fps0.mp4 -r 1 scene%04d.bmp  

### Export desired frames   
Frame by frame  
ffmpeg -r 1 -i /home/wanglab/Programs/Hourglass/build/Phox2b_#14_teLC_Choco_20220908_3_350fps1_labeled.mp4 -vf trim=start_frame=40000:end_frame=45000 -r 1 scene%04d.bmp  

Frame by frame + increase brightness  
fmpeg -r 1 -i /home/wanglab/Programs/Hourglass/build/Phox2B#8_20220810_3_1_labeled.mp4 -vf eq=brightness=0.3:contrast=1.5 -r 1 scene%04d.png  

Every 2 frames  
ffmpeg -i /home/wanglab/Desktop/build-CameraViewer-Qt_Static-Debug/Phox2b_#14_teLC_Sprinkles_20220919_1_500fps1.mp4 -vf "select=not(mod(n\,2))" -vsync vfr scene%04d.png

__bmp format produces a better image quality.__

### Create a movie from frames  
ffmpeg -r 35 -f image2 -s 640X480 -start_number 51300 -i scene%04d.bmp -vframes 4000 -vcodec libx264 -crf 10 -pix_fmt yuv420p test.mp4  
or  
ffmpeg -r 35 -f image2 -s 640X480 -i scene%04d.bmp -vcodec libx264 -crf 10 -pix_fmt yuv420p test.mp4  
See: https://hamelot.io/visualization/using-ffmpeg-to-convert-a-set-of-images-into-a-video/  

-r is the framerate (fps).  
-crf is the quality, lower means better quality, 15-25 is usually good.  
-s is the resolution.  
-pix_fmt yuv420p specifies the pixel format, change this as needed.  

# Concatenate videos
make a text file 'video.txt'
List the files to concatenate like this:
file 'Phox2b_#19_teLC_20221205_1_350fps1.mp4'
file 'Phox2b_#19_teLC_20221205_2_350fps1.mp4'

Then run this command from the directory that contains video files and video.txt:
ffmpeg -f concat -safe 0 -i video.txt -c copy output.mp4

### Slowing down an exsiting video  
ffmpeg -y -i /home/wanglab/Programs/Hourglass/build/Phox2B#8_20220720_1_1_labeled.mp4 -r 35 -s 640x480 -c:v libx264 -b:v 3M -strict -2 -movflags faststart /home/wanglab/Programs/Hourglass/build/Phox2B#8_20220720_1_1_labeled_slow.mp4   

 ## How to track the tongue and jaw
 1. Open the licking_jaw_config.json (/home/wanglab/Programs/Hourglass/build/licking_jaw_config.json)
 2. Specify the video file to be tracked.  
    "videos": "/media/wanglab/jun/Licking/Phox2b_videos/Phox2b#8_IRt_TeLC/Phox2b_#8_teLC_20220808_1_350fps1.mp4",   
 3. Open Terminal
 4. Change directory: /home/wanglab/Programs/Hourglass/build/  
 5. Run: ./Hourglass -d licking_jaw_config.json --train=false  
    Check -h help to see options
 6. Two .5 and one labeled video will be generated. 
  
## extract_h5.m 
**_For video tracking data_**  
Matlab function. Extract information from h5 files. 
Source file: h5 files generated by the tongue-jaw tracker tat Paul generated.  
Tongue: 0.h5   
Jaw: 1.h5  

## NumberOfLicksPerBout.m
**_For video tracking data_**  
For analyzing licking data from tracked video in Matlab.

## Density_plot_Python_ver2.ipynb
**_For data produced in cell_count_**  
For plotting scatter and density plots of neurons. It was originally written in Julia. It has been converted to Python.

## LickParameters_Python_20230322.ipynb
**_For lickometer data_**  
Python code:  This code analyze licking parameters from mat files generated by "Licking_Optogenetics_LEDCue.m" and save them in a csv file.
The parameters analyzed by this code are: "Licks/Bout", "MultiLicks/Bout", "Bout Duration", "Frequency"

## Convert_SessionData.m
**_For lickometer data_**  
This code convert SessionData files generated by the behavior padadigm "LickToGetReward_Traing" and "LickToGetReward_5" into
adjusted_SessionData files. In the adjusted files Cue start is adjusted to 0. This is for the processing in Python. Updated: 2023-04-18

## Plot_LickToGetReward_2023-04-19.ipynb & Plot_LickToGetReward_batch_2023-04-19.ipynb
**_For lickometer data_**  
Plot data from LickToGetReward. Before running these codes, you need to process SessionData by Convert_SessionData.m.

## extract_lick_parameters_from_LickToGetReward_Training_2023-05-13.ipynb
**_For lickometer data_**    
Extract lick_parameters from LickToGetReward_Training
This script extract licking parameters, including "mean_num_licksperbout", "mean_num_multi_licks", "mean_bout_duration", 
"mean_frequency", "Avg. Licks in Dringking period ", "LickThreshold". 
This script contains "Count_licks_in_Drinking.ipynb"
