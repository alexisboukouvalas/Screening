function [meanEE, meanStarEE, stdEE, EE] = Morris(R,BallTraj, funcEvalArray, delta, linearTransformation, PstarTraj)
%% function to evaluate Morris design - pass in
%% Morris design (BallTraj) and simulator function evaluated at these
%% points - should be (k+1) * R size..
%% Also pass in delta (in [0,1] space) and linear transformation
%% Returns moments of Elementary Effect distribution.
%% mean, mean of absolute values, standard deviation and raw elementary
%% effects.

k = size(BallTraj,2);
EE = nan(R,k); % elementary effects

% Reverse linear transformation
% not available design is in [0,1]
if(exist('linearTransformation','var'))
    if(size(linearTransformation,1) ~= k || size(linearTransformation,2) ~= 2)
        error('Morris::linearTransformation should be k X 2 matrix');
    end
    for iFactor = 1:k
        if(linearTransformation(iFactor,2) < linearTransformation(iFactor,1))
            error('Morris::Min bigger than max for factor %g',iFactor);
        end
        BallTraj(:,iFactor) = (BallTraj(:,iFactor) -  linearTransformation(iFactor,1) ) / ( linearTransformation(iFactor,2) - linearTransformation(iFactor,1) );
    end
else
    % check 0,1 assumption
    if(sum(sum(BallTraj < 0)) > 0 || sum(sum(BallTraj > 1)) > 0)
        error('Morris:: Not in [0,1] range');
    end
end

for r=1:R
    idxStart = (k)*(r-1)+1; idxEnd = idxStart + (k) - 1;
    Pstar = PstarTraj( idxStart:idxEnd, :);
    if(size(Pstar,1) ~= k); error('bad index'); end
    
    idxStartBB = (k+1)*(r-1)+1; 
    
    % calculate elementary effects for each variable
    for iFactor=2:(k+1)
        iEE = find(Pstar(iFactor-1,:) == 1); % permutation matrix tells us which EE is calculated
        if(length(iEE) ~= 1); error('exactly one'); end
        if(iEE < 1 || iEE > k); error('bad var'); end

        iTraj = idxStartBB + iFactor - 1;
        if(BallTraj(iTraj,iEE) - BallTraj(iTraj-1,iEE) > 0)
            % we add delta
            EE(r,iEE) = (funcEvalArray(iTraj) - funcEvalArray(iTraj-1)) / delta;
        elseif(BallTraj(iTraj,iEE) - BallTraj(iTraj-1,iEE) < 0)
            % we remove delta
            EE(r,iEE) = (funcEvalArray(iTraj-1) - funcEvalArray(iTraj)) / delta;
        else
            error('cannot be zero');
        end
    end
end

% And calculate moment of EE
meanEE = mean(EE); 
meanStarEE = mean(abs(EE));
stdEE = std(EE);


