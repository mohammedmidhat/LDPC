% mpdec
% Binary LDPC decoder using the Sum-Product Algorithm (SPA)
% With emulated fixed-point precsion

% The original implementation by : Alan Bao Jian ZHOU
% It was modified to :
% 1- Emulate fixed-point precsion
% 2- fix the bug when the LLR message --> check node = 0
% Modified : 9/28/2017

function [cHat, nIteration] = mpdec(H, rxLLR, nIterationMax, tanh_table, atanh_table)
    
    % Code structure
    % Caution: Move this part out for efficiency, as an independent preprocessing function.
    [nCheck, nBit] = size(H);
    nEdge = nnz(H);
    indexLinear = find(H);
    [indexCheck, indexBit] = find(H); % For sparse matrix, find() is fast.
    % The following initialization of the check nodes takes ONE THIRD of the running time!
    checkNode = struct('indexBitConnected', 0);
    check = repmat(checkNode, nCheck, 1);
    for iCheck = 1:nCheck
        check(iCheck).indexBitConnected = indexBit(indexCheck == iCheck); % Connections to check nodes
    end % for iCheck

    % Messages
    % Caution: Move the memory allocation of this part out for efficiency.
    messageChannel = fpt(rxLLR);
    messageBit = messageChannel; % Intialize bit messages
    msgb2ch = spalloc(nCheck, nBit, nEdge);
    msgch2b = spalloc(nCheck, nBit, nEdge);
    cHat = zeros(1,nBit);
    for nIteration = 1:nIterationMax
        msgb2ch(indexLinear) = fpt(fpt(messageBit(indexBit)).' - msgch2b(indexLinear));
        for iCheck = 1:nCheck
            for iBit = check(iCheck).indexBitConnected'
                product = fpt(prod(cell2mat(values(tanh_table,num2cell(fpt(full(msgb2ch(iCheck, setdiff(getfield(check(iCheck),'indexBitConnected'),iBit))/2)))))));
                msgch2b(iCheck, iBit) = fpt(2*cell2mat(values(atanh_table,{product})));
            end
        end
        messageBit = fpt(full(sum(msgch2b, 1) + messageChannel));
        % Hard decision, then exit or continue
        cHat(messageBit >= 0) = 0;
        cHat(messageBit < 0) = 1;
        if sum( mod(cHat * H', 2) ) == 0
            return;
        end
    end % for nIteration
    
end