function showVideo( M, frameSize, output )
%% Shows three videos: original, background and sparse parts


% original video
toplay = M;

no_frames = size(toplay,2);
width = frameSize(1);
length = frameSize(2);
video_ori = zeros(width,length,no_frames); % width*length*frames
for i=1:no_frames
    video_ori(:,:,i) = reshape(toplay(:,i),width,length);
end

% background term video
toplay = output.L;
no_frames = size(toplay,2);
width = frameSize(1);
length = frameSize(2);
video_L = zeros(width,length,no_frames); % width*length*frames
for i=1:no_frames
    video_L(:,:,i) = reshape(toplay(:,i),width,length);
end

% sparse term video
toplay = output.S;
no_frames = size(toplay, 2);
width = frameSize(1);
length = frameSize(2);
video_S = zeros(width,length,no_frames); % width*length*frames
for i=1:no_frames
    video_S(:,:,i) = reshape(toplay(:,i),width,length);
end


video_combine = [video_ori, video_L, abs(video_S)];
insert = zeros(width,10,no_frames)+ max(max(max(video_combine)));
video_combine = [video_ori, insert, video_L, insert, abs(video_S)];

figure();

% display video now 
for i=1:no_frames
    imagesc(video_combine(:,:,i)); 
    colormap('gray'); axis off;  % title('FW-T');
    pause(.01)
end

end

