function saveVideo( workingDir, frameDirName, outputFileName )
%% Creates and saves a video from given frame files
%
%   Author: Vahan Hovhannisyan, 2017.

imageNames = dir(fullfile(workingDir, frameDirName, '*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir, outputFileName));
%outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo);

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir, frameDirName, imageNames{ii}));
   writeVideo(outputVideo, img)
end

close(outputVideo);

end

