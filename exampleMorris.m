%% Generate morris design
%close all; 
clear all;
rand('state',42);
randn('state',42);

funcEval = 'polyfunForMorris2D'; % function to evaluate


% Produce design
%linearTransformation = [-100 0; 100 200]; % max min for each factor
linearTransformation = [0 1; 0 1];
deltaP = 1; p = 10; k = 2; R = 5;
[BallTraj, PstarTraj, delta] = MorrisDesign(R,k,p,deltaP, linearTransformation);

if(k==2)
    % can plot
        figDesign = figure; hold on; title('morris design');  xlabel('x_1'); ylabel('x_2'); 
        col = {'b','g','r','c','m','y','k'};
        B = BallTraj;
        for r=1:R
            % each trajectory is k+1
            Bstar = B( ((r-1) * (k+1) + 1) :  ((r-1) * (k+1) + k+1), : );
            if(size(Bstar,1) ~= (k+1)); erro('bad len'); end
        
            colT = col{mod(r,length(col))+1};
            plot(Bstar(:,1), Bstar(:,2),  sprintf('o-%s',colT), 'MarkerSize',10,'MarkerFaceColor',colT);
            text(Bstar(:,1), Bstar(:,2),  sprintf('%g',r), 'FontSize',20);
        end        
end


%% Evaluate simulator
funcEvalArray = nan(size(BallTraj,1),1); % (k+1)* R runs 
for iRun = 1:size(BallTraj,1)
    funcEvalArray(iRun) = feval(funcEval, BallTraj(iRun,:) );
end
    
%% Calculate Elementary Effects
[meanEE, meanStarEE, stdEE, EE] = MorrisEvaluate(R, BallTraj, funcEvalArray, delta, linearTransformation, PstarTraj);

    
% mean effect per variable
disp( [ meanEE' meanStarEE' stdEE'] );


figure; hold on;
xlabel('\mu^*'); ylabel('\sigma');
plot( meanStarEE, stdEE , 'ro', 'MarkerSize',10, 'MarkerFaceColor','r' );
for i=1:k
    text( meanStarEE(i), stdEE(i), sprintf('%g',i), 'FontSize', 30 );
end

