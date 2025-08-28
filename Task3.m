% 
% % Define parameters
% w = 5;% Size of the Gaussian filter for smoothing
% sigma = 1.5; % Standard deviation of the Gaussian filter for smoothing
% k = 0.04; % The constant used in the Harris response calculation
% p = 0.01; % The hyperparameter used to determine corner candidates
% 
% % Load the input image
% originalImage = imread('eight.tif');
% 
% % Ensure the image is grayscale
% if size(originalImage, 3) == 3
%     inputImage = rgb2gray(originalImage); % Convert RGB to grayscale
% 
% else
%     inputImage = originalImage;
% end
% 
% % convert image to double
% inputImage = im2double(inputImage);
% 
% % Call the Harris Corner Detector function
% ori_corners = harrisCornerDetector(inputImage, w, sigma, k, p);
% 
% 
% 
% intensityShiftedImage = imread('intensityShifted.png');
% 
% % Ensure the image is grayscale
% if size(intensityShiftedImage, 3) == 3
%     inputImage = rgb2gray(intensityShiftedImage); % Convert RGB to grayscale
% 
% else
%     inputImage = intensityShiftedImage;
% end
% 
% % convert image to double
% inputImage = im2double(inputImage);
% 
% intensityShifted_corners = harrisCornerDetector(inputImage, w, sigma, k, p);
% 
% 
% % plot the image with corners
% figure;
% 
% % Plot the original image 
% subplot(1, 2, 1);
% imshow(originalImage);         
% title('Original Image');
% hold on;
% [y, x] = find(ori_corners);          % Get row (y) and column (x) indices of the detected corners
% plot(x, y, 'ro', 'MarkerFaceColor', 'r');                % Plot the corners as red solid circle
% title('Original Image');
% hold off;
% 
% % Plot the corners overlaid on the image 
% subplot(1, 2, 2);
% imshow(intensityShiftedImage);              
% hold on;
% [y, x] = find(intensityShifted_corners);          % Get row (y) and column (x) indices of the detected corners
% plot(x, y, 'ro', 'MarkerFaceColor', 'r');                % Plot the corners as red solid circle
% title('Intensity Shifted Image');
% hold off;
% 
% 
% 
% 
% 
% function corners = harrisCornerDetector(inputImage, w, sigma, k, p)
% 
%     % Ensure filter size is odd
%     if mod(w, 2) == 0
%         error('Filter size w must be an odd integer.');
%     end
% 
%     %% Calculate derivatives
%     sobelX = [-1 0 1; -2 0 2; -1 0 1];
%     sobelY = [-1 -2 -1; 0 0 0; 1 2 1];
% 
%     Ix = imfilter(inputImage, sobelX, 'replicate');
%     Iy = imfilter(inputImage, sobelY, 'replicate');
% 
%     %% Calculate product of derivities
%     Ix2 = Ix.^2;
%     Iy2 = Iy.^2;
%     Ixy = Ix .* Iy;
% 
%     %% Apply Gaussian filter for smoothing
%     % generate gaussian filter with size w
%     gFilter = fspecial('gaussian', w, sigma);
% 
%     % Apply filter to derivities
%     A = imfilter(Ix2, gFilter, 'replicate');
%     B = imfilter(Ixy, gFilter, 'replicate');
%     C = imfilter(Iy2, gFilter, 'replicate');
% 
%     %% Calulate Harris Response, R
%     R = (A .* C - B.^2) - k * (A + C).^2;
% 
%     %% Extract corner candidate with above threshold
%     cornerCandidates = R > p * max(R(:)); % matrix of 0 and 1
% 
%     %% Extract local maxima
%     corners = false(size(R));
% 
%     % Loop through every row and column of cornerCandidates
%     for i = 2:size(cornerCandidates, 1)-1 % Avoid the borders
%         for j = 2:size(cornerCandidates, 2)-1
%             if cornerCandidates(i, j) % Check if the current pixel is a corner candidate
%                 % Extract the 3x3 neighborhood in R
%                 neighborhood = R(i-1:i+1, j-1:j+1);
% 
%                 % Check if the center pixel is the maximum in the neighborhood
%                 if R(i, j) == max(neighborhood(:)) && ...
%                    sum(neighborhood(:) == R(i, j)) == 1 % Ensure it's unique (no neighbor has the same R value, only the pixel itself)
%                     corners(i, j) = true; % Mark the pixel as a corner
%                 end
%             end
%         end
%     end
% end


% %% manipulate the image
% Read the original image
originalImage = imread('watermark_image.png');
% 
% % 1. Intensity Shift
% C = 50; % Constant value for shifting
% intensityShifted = originalImage + C; % Add C (or subtract for decreasing intensity)
% imwrite(intensityShifted, 'intensityShifted.png'); % Save the intensity-shifted image
% 
% % 2. Intensity Scaling
% k = 1.5; % Scaling factor
% intensityScaled = uint8(double(originalImage) * k); % Convert to double to avoid overflow
% imwrite(intensityScaled, 'intensityScaled.png'); % Save the intensity-scaled image
% 
% % 3. Image Translation
% dx = 20; % Translation in x-direction
% dy = 10; % Translation in y-direction
% translatedImage = imtranslate(originalImage, [dx, dy]);
% imwrite(translatedImage, 'translatedImage.png'); % Save the translated image
% 
% % 4. Spatial Scaling
s = 0.5; % Scaling factor (0.5 for downscaling)
spatiallyScaled = imresize(originalImage, s);
imwrite(spatiallyScaled, 'spatiallyScaled.png'); % Save the spatially scaled image
% 
% % Display confirmation
% disp('All manipulated images have been saved:');
% disp('1. intensityShifted.png');
% disp('2. intensityScaled.png');
% disp('3. translatedImage.png');
% disp('4. spatiallyScaled.png');
