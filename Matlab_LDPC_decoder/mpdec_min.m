% mpdec
% Binary LDPC decoder using the Min-Sum Algorithm (MSA)

% The original implementation by : Alan Bao Jian ZHOU
% It was modified to :
% 1- Implement the MSA instead of the SPA originally implemented
% Modified : 8/18/2017


function [cHat, nIteration, success] = mpdec_min(H, rxLLR, nIterationMax)
    
    success = 0;
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
    messageChannel = rxLLR;
    messageBit = messageChannel; % Intialize bit messages
    msgb2ch = spalloc(nCheck, nBit, nEdge);
    msgch2b = spalloc(nCheck, nBit, nEdge);
    cHat = zeros(1,nBit);
    for nIteration = 1:nIterationMax
        msgb2ch(indexLinear) = messageBit(indexBit).' - msgch2b(indexLinear);
        for iCheck = 1:nCheck
            for iBit = check(iCheck).indexBitConnected'
                msg = msgb2ch(iCheck, setdiff(getfield(check(iCheck),'indexBitConnected'),iBit));
                msgch2b(iCheck, iBit) = prod(sign(msg))*min(abs(msg));
            end
        end
        messageBit = sum(msgch2b, 1) + messageChannel;
        % Hard decision, then exit or continue
        cHat(messageBit >= 0) = 0;
        cHat(messageBit < 0) = 1;
        if sum( mod(cHat * H', 2) ) == 0
            success = 1;
            return;
        end
    end % for nIteration
    
end % fxn