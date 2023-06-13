% Get all h5 files in the current directory
files = dir('*.h5');

% Sort the files so that matching tongue and jaw files are adjacent
files = sort({files.name});

% Check if there are an even number of files (should be pairs of tongue and jaw files)
if mod(length(files),2) ~= 0
    error('The number of h5 files should be even (pairs of tongue and jaw files).');
end

% Process each pair of files
for i = 1:2:length(files)
    % Load tongue and jaw files
    Tongue = files{i};
    Jaw = files{i+1};

    % Continue with your existing code...
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
        frame_ind = j;
        
        %{
        heights{frame_ind}(heights{frame_ind} == 0) = 1;
        heights{frame_ind}(widths{frame_ind} == 0) = 1;
        
        img(heights{frame_ind},widths{frame_ind}) = probs{frame_ind};
        %}
        
        for k=1:length(heights{frame_ind})
            if heights{frame_ind}(k) == 0
                heights{frame_ind}(k) = 1;
            end
            if widths{frame_ind}(k) == 0
                widths{frame_ind}(k) = 1;
            end
            img(heights{frame_ind}(k),widths{frame_ind}(k)) = probs{frame_ind}(k);
        end
        

        img(img<0.2) = 0;

        areas(frames(j),1) = sum(sum(img));
        disp(['Processed frame: ', num2str(j)])
    end

    % Store areas and Jaw_heights
    S.areas = areas;
    S.Jaw_heights = Jaw_heights{1};

    % Save the S using the original (parent) file names
    [~,name,~] = fileparts(Tongue);
    name = name(1:end-2); % remove the last two characters (_0 or _1)
    save([name '.mat'], 'S');
end