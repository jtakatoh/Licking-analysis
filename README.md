# Licking-analysis

## install ffmpeg
https://cyber-tenchou.com/tips/20200904/  
https://eng-entrance.com/linux-vi-save  
vi .zshrc  
Write: export PATH=/Applications:$PATH  
To save and exit: :wq 

## Making masks

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
  
 ## Tongue_Jaw_correlation.m  
  Plot tongue area and jaw heigts.  
  Calculate correlation coefficient between tongue and jaw movements.  
  Source file: h5 files generated by the tongue-jaw tracker tat Paul generated.  
  Tongue: 0.h5   
  Jaw: 1.h5  
