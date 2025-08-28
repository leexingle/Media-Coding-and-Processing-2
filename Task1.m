%% call the function by extractHOGFeatures(inputImage) with parameter of an image 
% Load the input image
inputImage = imread('autumn.jpg');

% Ensure the image is grayscale
if size(inputImage, 3) == 3
    inputImage = rgb2gray(inputImage); % Convert RGB to grayscale
end

% Ensure image dimensions are divisible by 8, otherwise Resize the image
% Get size of image
[height, width] = size(inputImage);

% Compute the new dimensions that are divisible by 8
newHeight = ceil(height / 8) * 8;
newWidth = ceil(width / 8) * 8;

if height ~= newHeight || width ~= newWidth
    % Resize the image using imresize
    inputImage = imresize(inputImage, [newHeight, newWidth], 'bilinear');
end

% Normalize the image to double precision
inputImage = im2double(inputImage);

% Call the HOG feature extraction function
extractedHOGFeatures = extractHOGFeatures(inputImage);

% Display the length of the extracted features
disp(['Length of extracted HOG features: ',num2str(length(extractedHOGFeatures))]);


function extractedHOGFeatures = extractHOGFeatures (inputImage)
    %% calculate gradient and direction
    % Sobel filter for x and y direction
    Gx = [-1 0 1; -2 0 2; -1 0 1]; % Gradient in x-direction
    Gy = [-1 -2 -1; 0 0 0; 1 2 1]; % Gradient in y-direction

    % Compute gradients
    dx = filter2(Gx, inputImage, 'same'); % df/dx
    dy = filter2(Gy, inputImage, 'same'); % df/dy

    % calculate magnitude and direction
    gradientMagnitude = abs(dx) + abs(dy);
    gradientDirection = atan2(dy, dx); % in radian
    
    % Convert angle from radian to degrees 
    gradientDirection = rad2deg(gradientDirection);
    % angle normalize to [0, 180]
    gradientDirection(gradientDirection < 0) = gradientDirection(gradientDirection < 0) + 180;

    %% Divide image into 8x8 cells
    cellSize = 8; 
    [rows, cols] = size(inputImage);
    numCellsX = floor(cols / cellSize);
    numCellsY = floor(rows / cellSize);

    %% Calulate histogram of gradients
    numBins = 9;
    HOG = zeros(numCellsY, numCellsX, numBins); % matrix to store histograms cells

    % Calculate histogram for each cell
    for i = 1:numCellsY
        for j = 1:numCellsX

            % Define cell boundaries
            startRow = (i-1) * cellSize + 1;
            endRow = i * cellSize;
            startCol = (j-1) * cellSize + 1;
            endCol = j * cellSize;

            % Extract cell gradients
            cellMagnitude = gradientMagnitude(startRow:endRow, startCol:endCol);
            cellDirection = gradientDirection(startRow:endRow, startCol:endCol);

            % Compute histogram for the cell
            histValues = zeros(1, numBins);
            binEdges = linspace(0, 180, numBins + 1);
            
            % Loop through each pixel in the cell
            for i = 1:cellSize
                for j = 1:cellSize
                    % Get the gradient orientation for the current pixel
                    angle = cellDirection(i, j);
                    
                    % Find the two nearest bin edges
                    binIndex1 = find(angle >= binEdges(1:end-1) & angle < binEdges(2:end), 1);
                    binIndex2 = binIndex1 + 1; % The next bin
                    % Ensure binIndex2 does not exceed the number of bins
                    if binIndex2 > numBins
                        binIndex2 = 1;  % Wrap around to the first bin if necessary
                    end

                    % Check if angle is exactly on a bin edge
                    if angle == binEdges(binIndex1)
                        % If angle is exactly on a bin edge, assign the full magnitude to that bin
                        histValues(binIndex1) = histValues(binIndex1) + cellMagnitude(i, j);
                    else
                        % Calculate the distance to both bin edges
                        dist1 = abs(angle - binEdges(binIndex1));
                        dist2 = abs(angle - binEdges(binIndex2));
                        
                        % Normalize the distances to find the proportional contribution
                        totalDist = dist1 + dist2;
                        
                        % Assign proportion of the magnitude to each bin based on distance
                        histValues(binIndex1) = histValues(binIndex1) + (dist2 / totalDist) * cellMagnitude(i, j);
                        histValues(binIndex2) = histValues(binIndex2) + (dist1 / totalDist) * cellMagnitude(i, j);
                    end
                end
            end

            % Store the histogram for the cell in the HOG matrix
            HOG(i, j, :) = histValues;
        end
    end

    %% Block Normalization
    blockSize = [2, 2];  % Block size: 2x2 cells
    numBlocksY = numCellsY - blockSize(1) + 1;  % Number of blocks in the vertical direction
    numBlocksX = numCellsX - blockSize(2) + 1;  % Number of blocks in the horizontal direction
    
    % Preallocate matrix for normalized HOG descriptors
    normalizedHOG = []; % list to store the whole normalised HOG
    
    % Loop through each possible block position
    for i = 1:numBlocksY
        for j = 1:numBlocksX
            % Concatenate histograms of the 2x2 block of cells
            blockHist = []; % store concatenated histogram of a block
            
            for m = 0:blockSize(1)-1  % Loop through rows of the block
                for n = 0:blockSize(2)-1  % Loop through columns of the block
                    % Extract the histogram for the cell (i+m, j+n)
                    cellHist = squeeze(HOG(i+m, j+n, :));  % Get histogram of the current cell into an array
                    blockHist = [blockHist; cellHist];  % Concatenate cellHist to the block histogram
                end
            end
            
            % L2 Normalization of the block histogram (sqrt of sum of squares)
            normBlockHist = blockHist / sqrt(sum(blockHist.^2) + eps);  % eps to avoid division by zero
            
            % Append the normalized block histogram to the final HOG descriptor
            normalizedHOG = [normalizedHOG; normBlockHist];
        end
    end
   
    %% Form HOG feature vector
    extractedHOGFeatures = reshape(normalizedHOG, [], 1); % concatenate into one vector
end
