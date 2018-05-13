% channel_cap
% Compute the capacity of the flash channel as captured by its
% confusion matrix. Note this is the capacity for a specific coding
% associated with the matrix
% Arguments:
% Output:
% cap: capacity (<= 4 (for QLC))
% shannon_limit: normalized capacity (the argument name is not very dscriptive, sorry)

% Author: Mohammed Al Ai Baky
% Created: 9/24/2017

function [cap shannon_limit] = channel_cap(conf_mat)
    cap = 0;

    conf_mat_size = size(conf_mat);
    y_size = conf_mat_size(1);
    x_size = conf_mat_size(2);
    
    P_x = 1/x_size;
    P_y = zeros(1,y_size);
    for y = 1:y_size
        for x = 1:x_size
            P_y(1,y) = P_y(1,y) + conf_mat(y,x)*P_x;
        end
    end
    
    for y = 1:y_size
        for x = 1:x_size
            P_xy = conf_mat(y,x)*P_x;
            if(P_xy)
                cap = cap + P_xy*log2(P_xy/(P_x*P_y(y)));
            end
        end
    end
    
    shannon_limit = cap/log2(x_size);
    
    return;
end
