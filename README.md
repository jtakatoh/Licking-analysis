# Licking-analysis

## Making masks
Output frames from vides  
Check te number of total frames  
**ffprobe -v error -select_streams v:0 -count_packets -show_entries stream=nb_read_packets -of csv=p=0 /home/wanglab/Desktop/build-CameraViewer-Qt_Static-Debug/Phox2b_#14_teLC_Choco_20220908_1_350fps0.mp4**  

VLC can output frames but it drops frames. ffmpeg is better.  
Change directory to were you wwant to output frame images and run this  
**ffmpeg -r 1 -i /home/wanglab/Desktop/build-CameraViewer-Qt_Static-Debug/Phox2b_#14_teLC_Choco_20220908_3_350fps0.mp4 -r 1 scene%04d.png**  

Create a movie from frames  
**ffmpeg -r 35 -f image2 -s 640X480 -start_number 51300 -i scene%04d.png -vframes 4000 -vcodec libx264 -crf 10 -pix_fmt yuv420p test.mp4**  
See: https://hamelot.io/visualization/using-ffmpeg-to-convert-a-set-of-images-into-a-video/  

 
 ### How to track the tongue and jaw
 1. Open the licking_jaw_config.json (/home/wanglab/Programs/Hourglass/build/licking_jaw_config.json)
 2. Specify the video file to be tracked.  
    "videos": "/media/wanglab/jun/Licking/Phox2b_videos/Phox2b#8_IRt_TeLC/Phox2b_#8_teLC_20220808_1_350fps1.mp4",   
 3. Open Terminal
 4. Change directory: /home/wanglab/Programs/Hourglass/build/  
 5. Run: ./Hourglass -d licking_jaw_config.json --train=false
    Check -h help to see options
 6. Two .5 and one labeled video will be generated. 
  
 ### Tongue_Jaw_correlation.m  
  Plot tongue area and jaw heigts.  
  Calculate correlation coefficient between tongue and jaw movements.  
  Source file: h5 files generated by the tongue-jaw tracker tat Paul generated.  
  Tongue: 0.h5   
  Jaw: 1.5  
