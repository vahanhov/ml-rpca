function [ D, imDims ] = loadImageFrames( rootPath, trainingDatabaseName, maxImages )
%%  Loads videos from image files
%
%   Author: Vahan Hovhannisyan, 2017.


userDirectoryContents = list_image_files(fullfile(rootPath, trainingDatabaseName));

if isempty(userDirectoryContents)
    error(['No image files were found! Check your paths; there should be images in ' fullfile(rootPath, trainingDatabaseName)]);
end

% TODO: this must be autumated or set for each dictionary
imDims = [];

numImages = length(userDirectoryContents);
if maxImages > 0
    numImages = min(numImages, maxImages);
end

disp(['Loading images from ' rootPath '...']);

for fileIndex = 1 : numImages
    imageName = userDirectoryContents{fileIndex};
    % disp(['Using image file ' imageName '...']);

    imageFileName = fullfile( rootPath, trainingDatabaseName, imageName );
    [ ~, ~, ext ] = fileparts(imageFileName);
    if strcmp(ext, '.avi')
        [ D, imDims ] = loadVideo( imageFileName, maxImages );
    else
        image = im2double(imread(imageFileName));
        if size(image, 3) > 1
            image = rgb2gray(image);
        end

        if isempty(imDims)
            imDims = size(image);
            D = zeros( prod(imDims), numImages );
        else
            if any(imDims ~= size(image))
                warning(['skipping image ' imageFileName ' due to dimension inconsistance.']);
                continue;
            end
        end
        
        imcol = reshape(image, [prod(imDims) 1]); %im(:);
        D(:, fileIndex) = imcol ./ max(imcol); %norm(imcol);
    end
end

disp(['Loaded ' num2str(size(D, 2)) ' images from ' num2str(trainingDatabaseName) '.']);

return;

