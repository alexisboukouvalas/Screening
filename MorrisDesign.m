function [BallTraj, PstarTraj, delta01] = MorrisDesign(R,k,p,deltaP, linearTransformation)
%% function to generate morris design - pass in number of trajectories (R),
%% input dimension (k)
%% Also pass discretisation level (p) and multiplicative factor of delta/(p-1).
%% Pass in linearTransformation to linearly rescale design. It's an array
%% k X 2 size specifying min,max for each factor. If not specified, then
%% design is return in 0,1 domain.
%% Returns design (BallTraj), PstarTraj containing which factor changes and delta used for each factor (before lin
%% transformation!).

if(mod(p,2) ~= 0)
    error('MorrisDesign::p should be even');
end
if(exist('deltaP','var') && ~isnan(deltaP) )
    delta01 = deltaP / (p-1); % we actually pass in multiplicative factor over grid
else
    delta01 = p/(2*(p-1));    % default suggested in Morris paper
end

if(exist('linearTransformation','var'))
    % UNDONE Need delta per vabiable
    % To do
    % check 0 1
    if( sum(linearTransformation(:,2)) ~= k || sum(linearTransformation(:,1)) ~= 0)
        error('UNDONE can only have 0 1 bounds for now');
    end
%      delta = delta01 ./ (linearTransformation(:,2) - linearTransformation(:,1) ); % rescale delta
    delta = delta01;
else
    delta = delta01;
end

     
% 
% for pf = 1:(p-1)
%     disp( pf/(p-1) );
% end
xlevels = linspace(0,1,p);
    
m = k + 1; % number of points in each trajectory

EE = nan(R,k); % elementary effects

BallTraj = [];
PstarTraj = [];

for r = 1:R
    % for each trajectory
    
    % deterministic design
    B = ones(m,k) - triu(ones(m,k));
    
    % now apply randomisation
    
    % Dstar is diagonal matrix with +-1 with equal probability
    DstarDiag = ones(k,1);
    idx = randperm(k);
    DstarDiag (idx(1: round(length(k/2))) ) = -1; % if f(x-D) is allowed,

    Dstar = diag(DstarDiag);

    % choose base value
    xstart = nan(k,1); % col vector
    for i=1:k
        xstartgrid = xlevels ( xlevels <= (1-delta) );
        idx = randperm(length(xstartgrid));
        xstart(i) = xlevels( idx(1) );
        %xstart(i) = idx(1) ;
    end
    
    % permutation matrix
    Pstar = eye(k);
    idx = randperm(k);
    Pstar = Pstar(:,idx);
    
    % random orientation of B
    Bstar = (ones(m,1) * xstart' + (delta/2) * ( (2*B - ones(m,k)) *Dstar + ones(m,k))) * Pstar;
    
    % debug check design
    if(sum(sum(Bstar < 0)) > 0 || sum(sum(Bstar > 1)) > 0)
        error('MorrisDesign::bad design');
    end
    
    BallTraj = [BallTraj; Bstar]; % store all trajectories
    PstarTraj = [PstarTraj; Pstar];
end


if(exist('linearTransformation','var'))
    % apply linear transformation
    if(size(linearTransformation,1) ~= k || size(linearTransformation,2) ~= 2)
        error('MorrisDesign::linearTransformation should be k X 2 matrix');
    end
    for iFactor = 1:k
        if(linearTransformation(iFactor,2) < linearTransformation(iFactor,1))
            error('MorrisDesign::Min bigger than max for factor %g',iFactor);
        end
        BallTraj(:,iFactor) = BallTraj(:,iFactor) * ( linearTransformation(iFactor,2) - linearTransformation(iFactor,1) ) + linearTransformation(iFactor,1);
    end
    
end

% check matrices right size
if(size(BallTraj,1) ~= (k+1)*R )
    error('MorrisDesign::incorrect sizes');
end

