% fpt
% Converts floating point to fixed point
% Arguments:
% Input:
% x: double
% y: double, emulate fixed point

% Author: Mohammed Al Ai Baky
% Created: 9/28/2017

function y = fpt(x)
    
    y = double(sfi(x, 8, 4));

end % fxn