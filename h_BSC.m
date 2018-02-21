% h_BSC
% Compute the entropy of a BSC channel

% Author: Mohammed Al Ai Baky
% Created: 2/17/2018

function entropy = h_BSC(p)
    entropy = 1 + p*log2(p)+(1-p)*log2(1-p);
    
    return;
end
