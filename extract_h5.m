function [areas, Jaw_heights] = extract_h5()

% This function reads h5 files from the tongue jaw tracking and outputs
% tongue areas and jaw heights

    [file1,path1] = uigetfile('*_0.h5','Select the first .h5 file'); % the file ends with 0 is tongue 
    [file2,path2] = uigetfile('*_1.h5','Select the second .h5 file'); % the file ends with 1 is jaw
    
    if isequal(file1,0) || isequal(file2,0)
        error('User cancelled the file selection.');
    end
    
    Tongue = fullfile(path1,file1);
    Jaw = fullfile(path2,file2);
    
    image_height = 480;
    image_width = 640;

    frames=h5read(Tongue,'/frames');
    heights=h5read(Tongue,'/heights');
    widths=h5read(Tongue,'/widths');
    probs=h5read(Tongue,'/probs');

    Jaw_heights=h5read(Jaw,'/heights');
    Jaw_widths=h5read(Jaw,'/widths');

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

        %img_array(j) = {img}; %This stores entire image array and uses
        %lots of memory
        areas(frames(j),1) = sum(sum(img));
        disp(['Processed frame: ', num2str(j)])
    end
end