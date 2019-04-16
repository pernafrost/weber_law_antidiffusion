

% This script simulates the movement of particles across a stable
% distribution, such as e.g. a grayscale image
% As a single image does not in general present all the combinations of
% adjacent gray levels, here I run the simulation on multiple stable
% density landscapes and put together the results. Each density landscape
% is basically low-pass filtered two-dimensional noise.
% 
% Once the landscape is created, I want to explore movements of particles
% that leave the density landscape unchanged. I achieve this by first 
% calculating the flows of particles that would be obtained during one time
% step through diffusion alone, and then immediately compensating them with
% either the same particles coming back, or other particles from the
% destination moving in the opposite direction.
% For each node (pixel) there are four edges towards the four adjacent
% neighbours. Here I implement the movement of particles through the
% lattice in the following way
% 1) Diffusion step: I choose a node with a probability proportional to the
% number of individuals in the node
% 2) I pick one random individual from the node and make it move to any of
% the adjacent nodes (including those with zero density, that is, zero
% probability of moving there)
% 3) I now pick a random individual from the destination node (which could
% also be the same individual that I just moved) and move it to the origin
% node. If the destination node had zero density, then I have 100% chance
% of picking the same individual that I just moved and putting it back.
%
% With the three steps above, the flow from node i to j is F_iforward = D * N_i, but
% then I have to subtract the number of particles that I put back to i among those
% that I just moved from i, which is the product of how many particles I put back
% and of their fraction in the total population at destination:
% F_ibackward = F_iforward * F_iforward / (F_iforward + N_j)
% which gives F_ij_by_diffusion = D * N_i - (D * N_i)^2 / (D * N_i + N_j)
%
% The compensatory flow from j to i is equal to the flow by diffusion
% but in the opposite direction and involves particles that were in j moving to i
%
% F_ij = F_ij_by_diffusion + F_ij_compensatory
% and
% F_ji = F_ji_by_diffusion + F_ji_compensatory
%





%% Step 0: Cleaning and parameters

clear all; % clear variables
close all; % close figures
delete(findobj(allchild(0), '-regexp', 'Tag', '^Msgbox_')) % close message boxes

prefixToFileNames=''; % this is useful if we want to run the simulation again but give new names

% rng shuffle;
rng(700); % this is with the same randomization every time (for debugging and for reproducibility)

% Parameters of visualisation
useBlackBackground = 0; % this is just the background colour of figures

isSavingFigures = 1; 
imTitle = 'single_circle_diff'; % title of saved figures
imCounter = 1;
figureFolder = 'figures';
exportImageResolution = '-r300';

if isSavingFigures && ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end


if useBlackBackground == 1
    colordef('black');
    figureBackgroundColour = [0, 0, 0];
    markerFaceColour = [1, 0.5, 0.3];
    lineColour = [1, 1, 1];
    % currentColourMap = hot(1000);
    % currentColourMap = currentColourMap(end:-1:1, :);
    
else
    colordef('white');
    figureBackgroundColour = [1, 1, 1];
    markerFaceColour = [0.85, 0.85, 0.8];
    lineColour = [0, 0, 0];
    % currentColourMap = hot(1000);
end
useFontName = 'PT Sans';
useFontSize = 24;

imgMin = 0; % minimum accepted value for the image
imgMax = 2048; % maximum accepted value for the image

% parameters of landscape
landscapeSize = 256;


NiBins = linspace(imgMin, imgMax, 200);
NjBins = linspace(imgMin, imgMax, 200);

centerLocalDensityBins = NiBins(1:end-1) + diff(NiBins)/2;
centerTargetDensityBins = NjBins(1:end-1) + diff(NjBins)/2;

% initialize empty arrays
nCountsNiNjTotal = zeros(length(NiBins), length(NjBins));
nCountsNiNjTargetTotal = zeros(length(NiBins), length(NjBins));

simCounter = 0;
tic;

possibleSigmaValues = [2048, 1024, 512, 256, 128, 64, 32, 8, 4];

% matlabpool
for countSigma = 1:length(possibleSigmaValues) % this is the sigma of the filter in cycles per image side
    sigma = possibleSigmaValues(countSigma);
    for desiredDC = 2.^(0.2:0.2:(log(imgMax)/log(2)))
        for desiredRMSContrast = 2.^(1:1:(log(imgMax)/log(2)))
            
            fprintf(1, 'sigma = %d, desiredDC = %d, desiredRMSContrast = %d\n', sigma, desiredDC, desiredRMSContrast);
            
            % parameters of particles
            nVisibleParticles = 200;
            nSteps = 256; % 32768;
            colourMapForParticles = jet(nVisibleParticles);
            isPlottingTrajectories = 0;
            D = 0.02; % this is the diffusion coefficient < 0.25
            
            %% Step 1: Creates a 2-dimensional density landscape

            % initial random matrix
            im = rand(landscapeSize);
            
            
            % %%%%%%%%%%%%%%  create a Gaussian filter in Frequency %%%%%%%%%%%%%
            [X1, Y1] = meshgrid(-landscapeSize/2:landscapeSize/2-1);
            gFT = ifftshift(exp(-(((X1).^2+(Y1).^2)/(2*sigma*sigma))));
            % gFT(1,1) = 0; % removes DC
            gFT(end/2+1,:) = 0; % for symmetry
            gFT(:,end/2+1) = 0;
            
            
            
            % Fourier transform the original image (in this case the noise)
            DC = mean(im(:));
            % imNoDC = im - DC;
            % NFFT = 2.^nextpow2(size(im)); % Next power of 2 from length of y
            % imFT = fft2(imNoDC, max(NFFT), max(NFFT));
            % f = linspace(0,1,max(NFFT)/2+1);
            
            imFT =fft2(im);
            % figure, imagesc(abs(imFT)); colormap(gray);
            % figure, stem(1:49, log(abs(imFT(1,2:50))),'Marker', 'o', 'LineWidth', 2, 'Color', lineColour, 'MarkerFaceColor', markerFaceColour);
            % xlabel('Frequency (cycles / image)', 'FontName', useFontName, 'FontSize', useFontSize);
            % ylabel('log(amplitude)', 'FontName', useFontName, 'FontSize', useFontSize);
            
            
            
            
            % %%%%%%%%%%%%%%%%%%%%%%% Filter the image %%%%%%%%%%%%%%%%%%%%%%%%%
            imFiltered = real(ifft2(imFT .* gFT));
            % figure, imagesc(imFiltered); colormap(gray); colorbar;
            
            
            % Now set scale and DC
            currentDC = mean(mean(imFiltered))
            currentRMSContrast = sqrt(sum(sum((imFiltered - currentDC).^2))/(landscapeSize*landscapeSize));
            imFinal = (imFiltered - currentDC)*desiredRMSContrast/currentRMSContrast + desiredDC;
            
            % Impose that the image should stay within the limits
            imFinal(imFinal < imgMin) = imgMin;
            imFinal(imFinal > imgMax) = imgMax;
            
            
            % imFinal = (imFinal > 400) * 400.00; % this is a check that
            % everything is fine: if particles move to the black region it means that
            % there is something wrong
            
            if (min(min(imFinal))) < 0
                errordlg( 'try to change the desired DC and rms contrast', 'negative densities');
            end
            
            % figure;
            % imagesc(imFinal); colormap(gray);
            % colorbar;
            
            % imFinal = double(rgb2gray(imread('Dora.png')))+60;
            
            

            nodeLeft = circshift(imFinal, [0, 1]);
            nodeRight = circshift(imFinal, [0, -1]);
            nodeUp = circshift(imFinal, [1, 0]);
            nodeDown = circshift(imFinal, [-1, 0]);
            nodeThis = imFinal;
            
            flowLeftByDiffusion = D * nodeThis - (D * nodeThis).^2 ./ (D * nodeThis + nodeLeft);
            flowLeftByDiffusion(nodeThis == 0 & nodeLeft == 0) = 0;
            
            flowRightByDiffusion = D * nodeThis - (D * nodeThis).^2 ./ (D * nodeThis + nodeRight);
            flowRightByDiffusion(nodeThis == 0 & nodeRight == 0) = 0;
            
            flowUpByDiffusion = D * nodeThis - (D * nodeThis).^2 ./ (D * nodeThis + nodeUp);
            flowUpByDiffusion(nodeThis == 0 & nodeUp == 0) = 0;
            
            flowDownByDiffusion = D * nodeThis - (D * nodeThis).^2 ./ (D * nodeThis + nodeDown);
            flowDownByDiffusion(nodeThis == 0 & nodeDown == 0) = 0;
            
            
            flowRightCompensatory = circshift(flowLeftByDiffusion, [0, -1]);
            flowLeftCompensatory = circshift(flowRightByDiffusion, [0, 1]);
            flowDownCompensatory = circshift(flowUpByDiffusion, [-1, 0]);
            flowUpCompensatory = circshift(flowDownByDiffusion, [1, 0]);
            
            flowLeft = flowLeftByDiffusion + flowLeftCompensatory;
            flowRight = flowRightByDiffusion + flowRightCompensatory;
            flowUp = flowUpByDiffusion + flowUpCompensatory;
            flowDown = flowDownByDiffusion + flowDownCompensatory;
            
            % % Test that the map is steady with the selected flows
            % im2 = imFinal - flowLeft + circshift(flowLeft, [0, -1])...
            %     - flowRight + circshift(flowRight, [0, 1]) ...
            %     - flowUp + circshift(flowUp, [-1, 0]) ...
            %     - flowDown + circshift(flowDown, [1, 0]);
            %
            % max(max(abs(imFinal - im2)))
            
            
            % Convert flows to probabilities (a same flow can result in different
            % probabilities from the point of view of a particle, for instance if a
            % node with very few particles must support a large flow, each particle
            % will have a comparatively high probability of moving).
            pLeft = flowLeft ./ imFinal;
            pRight = flowRight ./ imFinal;
            pUp = flowUp ./ imFinal;
            pDown = flowDown ./ imFinal;
            % these probabilities can be larger than one because D was only diffusion,
            % but then particles are required to move in order to compensate for events
            % in neighbouring cells
            % Here I check that the probabilities are not too high
            
            maxProb = max(max(pLeft + pRight + pUp + pDown))
            if maxProb > 1 % I stop if the probabilities are too high
                myMessage = sprintf('please try a smaller value of diffusion coefficient D.\nCurrent value of maxProb = %0.2f.\nCurrentvalue of D = %0.2f.', maxProb, D);
                errordlg(myMessage, 'D coefficient too big for this image');
                continue
            else
                simCounter = simCounter + 1;
            end
            
            % The alternative option would be to rescale the probabilities
            % pLeft = pLeft / maxProb;
            % pRight = pRight / maxProb;
            % pUp = pUp / maxProb;
            % pDown = pDown / maxProb;
            
            
            
            % build a pMap with all the probabilities together:
            pMap = cat(3, cumsum(cat(3, pLeft, pRight, pUp, pDown), 3), ones(size(pUp)));
            % now pMap (x,y,:) should have monotonically increasing numbers
            
            figDensityLandscape = figure;
            imagesc(imFinal); colormap(gray);
            cb = colorbar;
            set(gca, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14);
            set(cb, 'LineWidth', 1.5, 'TickDir', 'out', 'FontSize', 14);
            hold on;
            
            
            
            %% Now simulate a random walk of some selected particles on the map
            
            selectedParticlePositions = nan(nVisibleParticles, 2);
            
            
            % choose starting positions (proportional to the local density of the map)
            imHeight = size(imFinal, 1); imWidth = size(imFinal, 2);
            cumsumGrayLevels = reshape(cumsum(imFinal(:)), imHeight, imWidth);
            
            if isfinite(nVisibleParticles) && ceil(cumsumGrayLevels(end)) >= 1
                selectedParticles = randi(ceil(cumsumGrayLevels(end)), [nVisibleParticles, 1]);
                for pp = 1:nVisibleParticles
                    [selectedParticlePositions(pp, 1),selectedParticlePositions(pp, 2)] = find(cumsumGrayLevels >= selectedParticles(pp), 1, 'first');
                end
            end
            
            if isPlottingTrajectories
                for pp = 1:nVisibleParticles
                    plot(selectedParticlePositions(pp,2), selectedParticlePositions(pp,1), '.', 'MarkerSize', 22, 'Color', colourMapForParticles(pp,:));
                end
            end
            
            % Now make the particles walk and record the trajectories
            trajectoriesX = nan(nVisibleParticles, nSteps+1);
            trajectoriesY = nan(nVisibleParticles, nSteps+1);
            localDensity = nan(nVisibleParticles, nSteps+1);
            
            densityUp = nan(nVisibleParticles, nSteps+1);
            densityDown = nan(nVisibleParticles, nSteps+1);
            densityLeft = nan(nVisibleParticles, nSteps+1);
            densityRight = nan(nVisibleParticles, nSteps+1);
            movingDirection = nan(nVisibleParticles, nSteps);
            
            targetDensity = nan(nVisibleParticles, nSteps + 1); % the last time step will always be NaN because we do not move the particle after its end position
            nonTargetDensity1 = nan(nVisibleParticles, nSteps + 1);
            nonTargetDensity2 = nan(nVisibleParticles, nSteps + 1);
            nonTargetDensity3 = nan(nVisibleParticles, nSteps + 1);
            nonTargetDensity4 = nan(nVisibleParticles, nSteps + 1);
            
            trajectoriesX(:,1) = selectedParticlePositions(:,2); % be careful because sub2ind returns first the y index
            trajectoriesY(:,1) = selectedParticlePositions(:,1);
            localDensity(:,1) = imFinal(sub2ind(size(imFinal), trajectoriesY(:,1), trajectoriesX(:,1)));
            
            densityUp(:,1) = imFinal(sub2ind(size(imFinal), mod(trajectoriesY(:, 1) -1 - 1, size(imFinal,1)) + 1, trajectoriesX(:,1)));
            densityDown(:,1) = imFinal(sub2ind(size(imFinal), mod(trajectoriesY(:, 1) -1 + 1, size(imFinal,1)) + 1, trajectoriesX(:,1)));
            densityLeft(:,1) = imFinal(sub2ind(size(imFinal), trajectoriesY(:,1), mod(trajectoriesX(:, 1) - 1 - 1, size(imFinal,2)) + 1));
            densityRight(:,1) = imFinal(sub2ind(size(imFinal), trajectoriesY(:,1), mod(trajectoriesX(:, 1) +1 - 1, size(imFinal,2)) + 1));
            
            
            % Build a map based on trajectories and check that it converges to the
            % original map. This is just for checking that everything is fine.
            trajCumulativeMap = zeros(size(imFinal));
            
            for ss = 2:nSteps + 1
                for pp = 1:nVisibleParticles
                    % move the particle and record its position
                    direction = find(rand(1) <= pMap(trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1), :), 1, 'first');
                    movingDirection(pp,ss - 1) = direction; % stores direction
                    switch direction
                        case 1 % because of the way pMap is built, this is left
                            trajectoriesX(pp,ss) = mod(trajectoriesX(pp, ss-1) - 1 - 1, size(imFinal,2)) + 1; % I subtract 1 twice, once for moving left and once for the modulus operation
                            trajectoriesY(pp,ss) = trajectoriesY(pp, ss-1);
                            targetDensity(pp,ss-1) = nodeLeft(sub2ind(size(nodeLeft), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity1(pp,ss-1) = nodeRight(sub2ind(size(nodeRight), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity2(pp,ss-1) = nodeUp(sub2ind(size(nodeUp), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity3(pp,ss-1) = nodeDown(sub2ind(size(nodeDown), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity4(pp,ss-1) = NaN;
                        case 2 % because of the way pMap is built, this is right
                            trajectoriesX(pp,ss) = mod(trajectoriesX(pp, ss-1) +1 - 1, size(imFinal,2)) + 1;
                            trajectoriesY(pp,ss) = trajectoriesY(pp, ss-1);
                            nonTargetDensity1(pp,ss-1) = nodeLeft(sub2ind(size(nodeLeft), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            targetDensity(pp,ss-1) = nodeRight(sub2ind(size(nodeRight), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity2(pp,ss-1) = nodeUp(sub2ind(size(nodeUp), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity3(pp,ss-1) = nodeDown(sub2ind(size(nodeDown), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity4(pp,ss-1) = NaN;
                        case 3% because of the way pMap is built, this is up
                            trajectoriesX(pp,ss) = trajectoriesX(pp, ss-1);
                            trajectoriesY(pp,ss) = mod(trajectoriesY(pp, ss-1) -1 - 1, size(imFinal,1)) + 1;
                            nonTargetDensity1(pp,ss-1) = nodeLeft(sub2ind(size(nodeLeft), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity2(pp,ss-1) = nodeRight(sub2ind(size(nodeRight), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            targetDensity(pp,ss-1) = nodeUp(sub2ind(size(nodeUp), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity3(pp,ss-1) = nodeDown(sub2ind(size(nodeDown), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity4(pp,ss-1) = NaN;
                        case 4% because of the way pMap is built, this is down
                            trajectoriesX(pp,ss) = trajectoriesX(pp, ss-1);
                            trajectoriesY(pp,ss) = mod(trajectoriesY(pp, ss-1) -1 + 1, size(imFinal,1)) + 1;
                            nonTargetDensity1(pp,ss-1) = nodeLeft(sub2ind(size(nodeLeft), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity2(pp,ss-1) = nodeRight(sub2ind(size(nodeRight), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity3(pp,ss-1) = nodeUp(sub2ind(size(nodeUp), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            targetDensity(pp,ss-1) = nodeDown(sub2ind(size(nodeDown), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity4(pp,ss-1) = NaN;
                        case 5% because of the way pMap is built, this is stay
                            trajectoriesX(pp,ss) = trajectoriesX(pp, ss-1);
                            trajectoriesY(pp,ss) = trajectoriesY(pp, ss-1);
                            nonTargetDensity1(pp,ss-1) = nodeLeft(sub2ind(size(nodeLeft), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity2(pp,ss-1) = nodeRight(sub2ind(size(nodeRight), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity3(pp,ss-1) = nodeUp(sub2ind(size(nodeUp), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            nonTargetDensity4(pp,ss-1) = nodeDown(sub2ind(size(nodeDown), trajectoriesY(pp,ss-1), trajectoriesX(pp,ss-1)));
                            targetDensity(pp,ss-1) = NaN;
                    end
                    
                    trajCumulativeMap(trajectoriesY(pp,ss), trajectoriesX(pp,ss)) = trajCumulativeMap(trajectoriesY(pp,ss), trajectoriesX(pp,ss)) + 1;
                    
                end
                % record image density at current position
                localDensity(:,ss) = imFinal(sub2ind(size(imFinal), trajectoriesY(:,ss), trajectoriesX(:,ss)));
                % also record density around the current position
                densityUp(:,ss) = imFinal(sub2ind(size(imFinal), mod(trajectoriesY(:, ss) -1 - 1, size(imFinal,1)) + 1, trajectoriesX(:,ss)));
                densityDown(:,ss) = imFinal(sub2ind(size(imFinal), mod(trajectoriesY(:, ss) -1 + 1, size(imFinal,1)) + 1, trajectoriesX(:,ss)));
                densityLeft(:,ss) = imFinal(sub2ind(size(imFinal), trajectoriesY(:,ss), mod(trajectoriesX(:, ss) - 1 - 1, size(imFinal,2)) + 1));
                densityRight(:,ss) = imFinal(sub2ind(size(imFinal), trajectoriesY(:,ss), mod(trajectoriesX(:, ss) +1 - 1, size(imFinal,2)) + 1));
                
                
                
                if isPlottingTrajectories && nVisibleParticles > 0
                    cla;
                    imagesc(imFinal); colormap(gray);
                    % for the purpose of plotting I put nan to values on the edge to
                    % avoid long lines
                    trajectoriesXTemp = trajectoriesX;
                    trajectoriesXTemp(find(trajectoriesX == 1)) = NaN;
                    trajectoriesYTemp = trajectoriesY;
                    trajectoriesYTemp(find(trajectoriesY == 1)) = NaN;
                    for pp = 1:nVisibleParticles
                        plot(trajectoriesXTemp(pp,1:ss), trajectoriesYTemp(pp,1:ss), '-', 'LineWidth', 2, 'Color', colourMapForParticles(pp,:));
                        plot(trajectoriesX(pp,ss), trajectoriesY(pp,ss), '.', 'MarkerSize', 22, 'Color', colourMapForParticles(pp,:));
                    end
                end
                
                drawnow;
                
%                 if mod(ss, 100) == 0
%                     fprintf('step %d\n', ss);
%                 end
            end
            
            
            % If the simulation is correct, I expect that each node is visited in
            % proportion to the density of particles that it contains (as well as in
            % proportion to how many nodes with that particular number of particles
            % exist: p(particle is in node with density N) = #(nodes with density N) * (density N)
            numNodesWithDensityN = histc(imFinal(:), imgMin:1:imgMax);
            probOfBeingInDensityN = histc(localDensity(:), imgMin:1:imgMax);
            
            % this curve should be more or less flat:
            figure, plot(probOfBeingInDensityN ./ numNodesWithDensityN ./ (imgMin:1:imgMax)', '.');
            
            figure, imagesc(trajCumulativeMap); colormap(gray); colorbar; title('Cumulative map of all trajectories')
            
            
            
            %% Now analyse the trajectories in relation to encounteredDensities
            % First I collect all the data in terms of localDensity, destinationDensity
            % and whether the particle moved to destination.
            
            % put together the four densities around the current position
            % neighbourDensities = cat(3, densityLeft, densityRight, densityUp, densityDown);
            
            
            
            % The line below can be used to check that the vectors localDensity,
            % targetDensity etc. are correct
            % pp = 34; % choose one particle
            % disp([trajectoriesX(pp,1:end-1)', trajectoriesY(pp,1:end-1)', localDensity(pp,1:end-1)', densityLeft(pp,1:end-1)', densityRight(pp,1:end-1)', densityUp(pp,1:end-1)', densityDown(pp,1:end-1)', movingDirection(pp,:)', targetDensity(pp,:)', nonTargetDensity1(pp,:)', nonTargetDensity2(pp,:)', nonTargetDensity3(pp,:)', nonTargetDensity4(pp,:)'])
            % fprintf('X, Y, localDensity, left, right, up, down, direction, targetD, nonTargetD1, nonTargetD2, nonTargetD3, nonTargetD4\n')
            
            
            
            
            % bin the data putting together similar values of Ni and Nj (Ni = density
            % at the local position; Nj = density at the target destination.
            
            NjMinusNiBins = linspace(imgMin - imgMax, imgMax - imgMin, 201);
            NjMinusNiOverNiBins = linspace(-1, 1, 201);
            

            
            
            [nInNi, valInNi] = histc(localDensity, NiBins);
            [nInNjTarget, valInNjTarget] = histc(targetDensity, NjBins);
            [nInNjNonTarget1, valInNjNonTarget1] = histc(nonTargetDensity1, NjBins);
            [nInNjNonTarget2, valInNjNonTarget2] = histc(nonTargetDensity2, NjBins);
            [nInNjNonTarget3, valInNjNonTarget3] = histc(nonTargetDensity3, NjBins);
            [nInNjNonTarget4, valInNjNonTarget4] = histc(nonTargetDensity4, NjBins);
            
            [nInNjminusNiTarget, valInNjminusNiTarget] = histc(targetDensity - localDensity, NjMinusNiBins);
            [nInNjminusNiNonTarget1, valInNjminusNiNonTarget1] = histc(nonTargetDensity1 - localDensity, NjMinusNiBins);
            [nInNjminusNiNonTarget2, valInNjminusNiNonTarget2] = histc(nonTargetDensity2 - localDensity, NjMinusNiBins);
            [nInNjminusNiNonTarget3, valInNjminusNiNonTarget3] = histc(nonTargetDensity3 - localDensity, NjMinusNiBins);
            [nInNjminusNiNonTarget4, valInNjminusNiNonTarget4] = histc(nonTargetDensity4 - localDensity, NjMinusNiBins);
            
            
            clear nCountsNiNj nCountsNiNjTarget nCountsNiNjNonTarget freqMovedNiNj;
            
            % initialize empty arrays
            nCountsNiNj = zeros(length(NiBins), length(NjBins));
            nCountsNiNjTarget = zeros(length(NiBins), length(NjBins));
            nCountsNiNjNonTarget = zeros(length(NiBins), length(NjBins));
            freqMovedNiNj = NaN(length(NiBins), length(NjBins));
            
            
            % first construct the maps in DeltaX and DeltaY
            for cNi=min(valInNi(:)):max(valInNi(:)) % this should speed up a bit, because we 
                % cNi
                
                
                for cNj=1:length(NjBins) -1
                    
                    
                    % allNiVals = [];
                    % allNjVals = [];
                    allNCounts = 0;
                    allNCountsNjTarget = 0;
                    allNCountsNjNonTarget = 0;
                    
                    % Target
                    valNiNjT = intersect(find(valInNi==cNi), find(valInNjTarget==cNj));
                    
                    allNCounts = allNCounts + length(valNiNjT);
                    allNCountsNjTarget = allNCountsNjTarget + length(valNiNjT);
                    % allNCountsNjNonTarget = allNCountsNjNonTarget + 0;
                    % allNiVals = [allNiVals; localDensity(valNiNjT)];
                    % allNjVals = [allNjVals; targetDensity(valNiNjT)];
                    
                    % non target 1
                    valNiNjNT1 = intersect(find(valInNi==cNi), find(valInNjNonTarget1==cNj));
                    
                    allNCounts = allNCounts + length(valNiNjNT1);
                    % allNCountsNjTarget = allNCountsNjTarget + 0;
                    allNCountsNjNonTarget = allNCountsNjNonTarget + length(valNiNjNT1);
                    % allNiVals = [allNiVals; localDensity(valNiNjNT1)];
                    % allNjVals = [allNjVals; nonTargetDensity1(valNiNjNT1)];
                    
                    % non target 2
                    valNiNjNT2 = intersect(find(valInNi==cNi), find(valInNjNonTarget2==cNj));
                    
                    allNCounts = allNCounts + length(valNiNjNT2);
                    % allNCountsNjTarget = allNCountsNjTarget + 0;
                    allNCountsNjNonTarget = allNCountsNjNonTarget + length(valNiNjNT2);
                    % allNiVals = [allNiVals; localDensity(valNiNjNT2)];
                    % allNjVals = [allNjVals; nonTargetDensity1(valNiNjNT2)];
                    
                    % non target 3
                    valNiNjNT3 = intersect(find(valInNi==cNi), find(valInNjNonTarget3==cNj));
                    
                    allNCounts = allNCounts + length(valNiNjNT3);
                    % allNCountsNjTarget = allNCountsNjTarget + 0;
                    allNCountsNjNonTarget = allNCountsNjNonTarget + length(valNiNjNT3);
                    % allNiVals = [allNiVals; localDensity(valNiNjNT3)];
                    % allNjVals = [allNjVals; nonTargetDensity1(valNiNjNT3)];
                    
                    % non target 4
                    valNiNjNT4 = intersect(find(valInNi==cNi), find(valInNjNonTarget4==cNj));
                    
                    allNCounts = allNCounts + length(valNiNjNT4);
                    % allNCountsNjTarget = allNCountsNjTarget + 0;
                    allNCountsNjNonTarget = allNCountsNjNonTarget + length(valNiNjNT4);
                    % allNiVals = [allNiVals; localDensity(valNiNjNT4)];
                    % allNjVals = [allNjVals; nonTargetDensity1(valNiNjNT4)];
                    
                    
                    freqMovedNiNj(cNj,cNi) =  allNCountsNjTarget / allNCounts;
                    nCountsNiNj(cNj,cNi) =  allNCounts;
                    nCountsNiNjTarget(cNj,cNi) =  allNCountsNjTarget;
                    nCountsNiNjNonTarget(cNj,cNi) =  allNCountsNjNonTarget;
                    nCountsNiNjNonTarget(cNj,cNi) =  allNCountsNjNonTarget;
                end
            end
            
            
            nCountsNiNjTotal = nCountsNiNjTotal + nCountsNiNj;
            nCountsNiNjTargetTotal = nCountsNiNjTargetTotal + nCountsNiNjTarget;
            
            delete(findobj(allchild(0), '-regexp', 'Tag', '^Msgbox_'))
            close all;
            
            save('sim3_multiple_stable_landscapes_simulation.mat', 'nCountsNiNjTotal', 'nCountsNiNjTargetTotal', 'sigma', 'desiredDC', 'desiredRMSContrast', 'simCounter');
            fprintf(1, 'simCounter=%d, time = %f\n', simCounter, toc);
            
            % pause(60) % Only to let the computer cool down a bit
            
        end
    end
end




% create the colormap
x=(0:1:255)/255;

R = x/0.32 - 0.78125;
R(R<0) = 0;
R(R>1) = 1;

G = 2*x - 0.84;
G(G<0) = 0;
G(G>1) = 1;


B(find(x<=0.25)) = 4*x(x<=0.25);
B(intersect(find(x>0.25), find(x<=0.42))) = 1;
B(intersect(find(x>0.42), find(x<=0.92))) = -2*x(x>0.42 & x<=0.92) + 1.84;
B(find(x>0.92)) = x(x>0.92)/0.08 - 11.5;

myColourMap=[R', G', B'];




%% figure
cAxisLims = [0.0, 0.1];
figure,
set(gcf, 'Position', [1 1 1000 700]);
% imagesc(freqMovedNiNj);
h=pcolor(NiBins,NjBins,nCountsNiNjTargetTotal ./ nCountsNiNjTotal);
caxis(cAxisLims);
hold on;
axis equal xy tight;
box off;
colormap(myColourMap);
set(h,'edgecolor','none');
shading flat;

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');


lx = xlabel('N_i'); ly = ylabel('N_j');
set([lx ly], 'interpreter', 'tex');

colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
cly = ylabel(colorbarHandle, sprintf('P_{ij}'), 'FontSize', useFontSize);
set(cly, 'interpreter', 'tex');

colormap(myColourMap);

if isSavingFigures
    figureFileName = fullfile(figureFolder, [prefixToFileNames, 'probVsNiNj.eps']);
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end





%% figure slope in p vs. Ni

% make bins in Ni
% compute slope


figure,
set(gcf, 'Position', [1 1 1000 700]);
for column = 1:length(centerLocalDensityBins);

% slope in Pij = a (Nj - Ni), where a is proportional to 1/Ni
% now try to solve the system y = a x
% a = (Nj - Ni) \ pij
y = squeeze(((nCountsNiNjTargetTotal(1:(end -1), column))./ nCountsNiNjTotal(1:(end - 1), column)));
x = squeeze(centerTargetDensityBins - centerLocalDensityBins(column))';

notNaNValues = ~isnan(y);

slopeInPi1(column) = x(notNaNValues) \ y(notNaNValues)
linearFit(column,1:2) = polyfit(x(notNaNValues), y(notNaNValues), 1)
slopeInPi2(column) = (x(notNaNValues) - mean(x(notNaNValues))) \ (y(notNaNValues) - mean(y(notNaNValues)))

plot(x(notNaNValues) - mean(x(notNaNValues)), y(notNaNValues) - mean(y(notNaNValues)), '.');
hold on;
plot(x(notNaNValues) - mean(x(notNaNValues)), slopeInPi2(column) * (x(notNaNValues) - mean(x(notNaNValues))), 'r-');
title(['N_i in range ', num2str(round(NiBins(column))), ' - ', num2str(round(NiBins(column+1)))])
xlabel('N_j - N_i');
ylabel('p_{ij}')
pause(0.1)

set(gca, 'LineWidth', 1.5, 'box', 'off', 'FontName', 'Arial', 'FontSize', 20)

if isSavingFigures
    figureFileName = fullfile(figureFolder, [prefixToFileNames, 'pij_vs_NiNj', num2str(round(NiBins(column))), '-', num2str(round(NiBins(column+1))) '.eps']);
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end


cla
end




figure
% loglog(centerLocalDensityBins, slopeInPi2(:) .* centerLocalDensityBins', '.'); % this is constant
loglog(centerLocalDensityBins, slopeInPi2(:), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'Color', 'k', 'LineWidth', 1);
xLimits = get(gca, 'XLim');
yLimits = get(gca, 'YLim');
logScale = diff(yLimits)/diff(xLimits);  %# Scale between the x and y ranges
powerScale = diff(log10(yLimits))/...    %# Scale between the x and y powers
             diff(log10(xLimits));
set(gca,'DataAspectRatio', [1, logScale/powerScale, 1])
set(gca, 'LineWidth', 1.5, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');
lx = xlabel('N_i');
ly = ylabel('slope in P_{ij}', 'interpreter', 'tex');
set([lx ly], 'interpreter', 'tex');

hold on;
plot(centerLocalDensityBins, D ./ centerLocalDensityBins, 'k-', 'LineWidth', 2)


if isSavingFigures
    figureFileName = fullfile(figureFolder, [prefixToFileNames, 'slopePvsNi.eps']);
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end








%% figure
cAxisLims = [0.0, 0.1];
figure,
set(gcf, 'Position', [1 1 1000 700]);
% imagesc(freqMovedNiNj);
h=pcolor(NiBins,NjBins,nCountsNiNjTargetTotal ./ nCountsNiNjTotal);
caxis(cAxisLims);
hold on;
axis equal xy tight;
box off;
colormap(myColourMap);
set(h,'edgecolor','none');
shading flat;



for ii = 200:200:max(NiBins)
    plot(ii * ones(size(NjBins(5:5:end))), NjBins(5:5:end), 'k*');
end


set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');


lx = xlabel('N_i'); ly = ylabel('N_j');
set([lx ly], 'interpreter', 'tex');

colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
cly = ylabel(colorbarHandle, sprintf('P_{ij}'), 'FontSize', useFontSize);
set(cly, 'interpreter', 'tex');
colormap(myColourMap);

if isSavingFigures
    figureFileName = fullfile(figureFolder, [prefixToFileNames, 'probVsNiNjwithLines.eps']);
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end








%
% if isSavingFigures
%     figureFileName = sprintf('%s%cillustrative_drawing.eps', figureFolder, filesep);
%     set(gcf, 'InvertHardCopy', 'off');
%     set(gcf, 'PaperPositionMode', 'auto');
%     print(gcf, '-depsc2', '-tiff', '-r300', figureFileName); print(gcf, '-dpng', '-r300', [figureFileName(1:end-3), 'png'])
% end
%