function [imGray, pLeft, pRight, pUp, pDown] = perform_one_deblur_step(imGray, gamma, varargin)


if size(varargin, 2) >= 1
    isDeBlurringWithLaplacian = varargin{1};
else
    isDeBlurringWithLaplacian = 0;
end

if size(varargin, 2) >= 2
    meanField = varargin{2};
else
    meanField = 1;
end

if size(varargin, 2) >= 3
    avoidSaturationWhenDeblurring = varargin{3};
else
    avoidSaturationWhenDeblurring = 0;
end


if size(varargin, 2) >= 4
    iMax = varargin{4};
else
    iMax = inf;
end


if isDeBlurringWithLaplacian
    
    % I still need to compute the weber contrast for returning it and
    % for moving the visible particles
    
    
    % The max(0, value) indicates that particles only move up the gradient
    weberContrastRight = max(0, (imGray(:, [2:end, 1]) - imGray)./imGray);
    weberContrastLeft = max(0, (imGray(:, [end, 1:end - 1]) - imGray)./imGray);
    weberContrastDown = max(0, (imGray([2:end, 1], :) - imGray)./imGray);
    weberContrastUp = max(0, (imGray([end, 1:end - 1], :) - imGray)./imGray);
    
    % If there are no particles at some location, nothing moves out of
    % that location
    weberContrastLeft(imGray == 0) = 0;
    weberContrastRight(imGray == 0) = 0;
    weberContrastUp(imGray == 0) = 0;
    weberContrastDown(imGray == 0) = 0;
    
    %         weberContrastLeftRight = weberContrastLeft + weberContrastRight;
    %        weberContrastUpDown = weberContrastUp + weberContrastDown;
    %        weberContrastTot = weberContrastLeftRight + weberContrastUpDown;
    
    
    % the probability of NOT moving is given by the combined
    % probability of not moving in any direction
    pLeft = min(0.25, weberContrastLeft * gamma);
    pRight = min(0.25, weberContrastRight * gamma);
    pUp = min(0.25, weberContrastUp * gamma);
    pDown = min(0.25, weberContrastDown * gamma);
    
    
    
    % the laplacian kernel
    L = [0, 1, 0;...
        1,-4, 1;...
        0, 1, 0];
    imGray = wextend('2D','ppd',imGray,[1, 1]);
    if avoidSaturationWhenDeblurring
        
        % with the line below sum(imGray(:)) is not conserved because the
        % laplacian is computed on imGray .* (iMax - imGray)
    imGray = imGray - imGray .* (iMax - imGray) .* convn(imGray, L, 'same')*gamma /(iMax)^2;

    else
    imGray = imGray - convn(imGray, L, 'same')*gamma;
    end
    imGray = imGray(2:end-1, 2:end-1);
    
    
    
    
    
else
    
    
    
    % The max(0, value) indicates that particles only move up the gradient
    weberContrastRight = max(0, (imGray(:, [2:end, 1]) - imGray)./imGray);
    weberContrastLeft = max(0, (imGray(:, [end, 1:end - 1]) - imGray)./imGray);
    weberContrastDown = max(0, (imGray([2:end, 1], :) - imGray)./imGray);
    weberContrastUp = max(0, (imGray([end, 1:end - 1], :) - imGray)./imGray);
    
    % If there are no particles at some location, nothing moves out of
    % that location
    weberContrastLeft(imGray <= 0) = 0;
    weberContrastRight(imGray <= 0) = 0;
    weberContrastUp(imGray <= 0) = 0;
    weberContrastDown(imGray <= 0) = 0;
    
    %         weberContrastLeftRight = weberContrastLeft + weberContrastRight;
    %        weberContrastUpDown = weberContrastUp + weberContrastDown;
    %        weberContrastTot = weberContrastLeftRight + weberContrastUpDown;
    
    
    % the probability of NOT moving is given by the combined
    % probability of not moving in any direction
    pLeft = min(0.25, weberContrastLeft * gamma);
    pRight = min(0.25, weberContrastRight * gamma);
    pUp = min(0.25, weberContrastUp * gamma);
    pDown = min(0.25, weberContrastDown * gamma);
    
    
    if meanField
        % Here I move exactly the fraction of
        % particles that correspond to the probability
        
        % movingLeftUp = imGray .* pLeft .* pUp;
        % movingLeftDown = imGray .* pLeft .* pDown;
        if avoidSaturationWhenDeblurring
            movingLeft = min(imGray .* pLeft, (iMax - imGray(:, [end, 1:end - 1]))/4);
            movingRight = min(imGray .* pRight, (iMax - imGray(:, [2:end, 1]))/4);
            movingUp = min(imGray .* pUp, (iMax - imGray([end, 1:end - 1], :))/4);
            movingDown = min(imGray .* pDown, (iMax - imGray([2:end, 1], :))/4);
        else
            movingLeft = imGray .* pLeft;
            % movingRightUp = imGray .* pRight .* pUp;
            % movingRightDown = imGray .* pRight .* pDown;
            movingRight = imGray .* pRight;
            movingUp = imGray .* pUp;
            movingDown = imGray .* pDown;
        end
        
        totalMovingParticles = movingLeft + movingRight + movingUp + movingDown;% + ...
        %movingLeftUp + movingLeftDown + movingRightUp + movingRightDown;
        
    else
        
        % here I draw the number of particles moving in each direction
        imGrayRound = round(imGray);
        %
        
        if avoidSaturationWhenDeblurring
            movingLeft = min(binornd(imGrayRound, pLeft), (iMax - imGray(:, [end, 1:end - 1]))/4);
            % movingRightUp = imGray .* pRight .* pUp;
            % movingRightDown = imGray .* pRight .* pDown;
            movingRight = min(binornd(imGrayRound, pRight), (iMax - imGray(:, [2:end, 1]))/4);
            movingUp = min(binornd(imGrayRound, pUp), (iMax - imGray([end, 1:end - 1], :))/4);
            movingDown = min(binornd(imGrayRound, pDown), (iMax - imGray([2:end, 1], :))/4);
        else
            %         movingLeftUp = binornd(imGrayRound, pLeft .* pUp);
            %         movingLeftDown = binornd(imGrayRound, pLeft .* pDown);
            movingLeft = binornd(imGrayRound, pLeft);
            %         movingRightUp = binornd(imGrayRound, pRight .* pUp);
            %         movingRightDown = binornd(imGrayRound, pRight .* pDown);
            movingRight = binornd(imGrayRound, pRight);
            movingUp = binornd(imGrayRound, pUp);
            movingDown = binornd(imGrayRound, pDown);
        end
        totalMovingParticles = movingLeft + movingRight + movingUp + movingDown;% + ...
        %             movingLeftUp + movingLeftDown + movingRightUp + movingRightDown;
        %         % the problem is that somewhere there will be more particles that want to
        %         % move than the particles available. If this is the case, I rescale
        %         % everything
        tooManyParticles = find(totalMovingParticles > imGray);
        movingLeft(tooManyParticles) = (movingLeft(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        movingRight(tooManyParticles) = (movingRight(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        movingUp(tooManyParticles) = (movingUp(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        movingDown(tooManyParticles) = (movingDown(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        %         movingLeftUp(tooManyParticles) = floor(movingLeftUp(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        %         movingRightUp(tooManyParticles) = floor(movingRightUp(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        %         movingLeftDown(tooManyParticles) = floor(movingLeftDown(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        %         movingRightDown(tooManyParticles) = floor(movingRightDown(tooManyParticles) .* imGray(tooManyParticles) ./ totalMovingParticles(tooManyParticles));
        %
        %         % recompute the total number of moving particles, which has changed
        totalMovingParticles = movingLeft + movingRight + movingUp + movingDown;%  + ...
        %             movingLeftUp + movingLeftDown + movingRightUp + movingRightDown;
        
    end
    
    
    
    
    
    
    
    
    % check that the sum of moving particles in each direction is
    % always equal to the total number of moving particles (that the sum
    % below is zero)
    % sum(sum((movingLeft + movingRight + movingUp + movingDown - movingParticles).^2))
    
    imGray = max(0, imGray -totalMovingParticles + movingLeft(:, [2:end, 1]) + movingRight(:, [end, 1:end-1]) ...
        + movingUp([2:end, 1], :) + movingDown([end, 1:end-1], :));% I impose that imGray is non negative (approximations)
    %+ movingLeftUp([2:end, 1], [2:end, 1]) + movingLeftDown([end, 1:end-1], [2:end, 1]) ...
    %+ movingRightUp([2:end, 1], [end, 1:end-1]) + movingRightDown([end, 1:end-1], [end, 1:end-1]);
    
end

