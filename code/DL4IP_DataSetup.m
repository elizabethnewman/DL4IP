function[X] = DL4IP_DataSetup(options)

% load data to reconstruct
switch options.data.name
    case 'moon'
        img = imread('moon.tif');
        img = prepareImage(img,options);
        X   = double(img);

    case 'GirlWithPearlEarring'
        % https://en.wikipedia.org/wiki/Girl_with_a_Pearl_Earring#/media/File:1665_Girl_with_a_Pearl_Earring.jpg
        img = imread('GirlWithPearlEarring.jpeg');
        img = prepareImage(img,options);
        X   = double(img);

    case 'peppers'
        img = imread('peppers.png');
        img = prepareImage(img,options);
        X   = double(img);

    case 'tomography'
        % requires AIRToolsII
        X = phantomgallery('shepplogan',options.data.size(1));
end

end


%% helper functions
function[img] = prepareImage(img,options)
    if ~isempty(options.data.size)
        img = imresize(img,options.data.size);
    end

    if ~options.data.color && size(img,3) == 3
        img = rgb2gray(img);
    end
end

