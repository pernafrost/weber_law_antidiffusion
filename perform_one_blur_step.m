function [imGray] = perform_one_blur_step(imGray, D, varargin)

if size(varargin, 2) >= 1
    isBlurringWithLaplacian = varargin{1};
else
    isBlurringWithLaplacian = 0;
end



% The difference between blurring through the convolution with a Laplacian
% and blurring by simulating the movement of particles are (1) speed of
% computation, (2) the laplacian is a 'mean field', (3) moving particles
% ensures that the pixels have integer values, the laplacian does not
if isBlurringWithLaplacian
    % the laplacian kernel
    L = [0, 1, 0;...
    1,-4, 1;...
    0, 1, 0];
    % blur the image through convolution with Laplacian
    % there are several possible decisions that can be taken for
    % boundary conditions, but I think we don't care too much
    imGray = wextend('2D','ppd',imGray,[1, 1]);
    imGray = imGray + convn(imGray, L, 'same')*D;
    imGray = imGray(2:end-1, 2:end-1);
else
    % blur the image by simulating movements of particles
    movingParticles = binornd(imGray, D*4);
    movingLeftRight = binornd(movingParticles, 0.5); % around half of the particles move to the left or to the right
    movingUpDown = movingParticles - movingLeftRight; % the remaining half moves up or down
    
    movingLeft = binornd(movingLeftRight, 0.5); % around half of the particles moving left or right move left
    movingRight = movingLeftRight - movingLeft; % and the remaining ones move to the right
    
    movingUp = binornd(movingUpDown, 0.5);
    movingDown = movingUpDown - movingUp;
    
    if sum(isnan(movingParticles(:))) > 0
        disp('found nan in the number of moving particles; this can happen because the binornd has found a non integer number?');
    end
    
    imGray = imGray -movingParticles + movingLeft(:, [2:end, 1]) + movingRight(:, [end, 1:end-1]) ...
        + movingUp([2:end, 1], :) + movingDown([end, 1:end-1], :);
    
end


