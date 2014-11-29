%%%%% parameters %%%%%%%%
timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
lambda=0.5; % coefficient of the weighted length term L(phi)
strong_alpha = 0.75;  % coefficient of the weighted area term A(phi)
alfa=strong_alpha; 
epsilon=1; % papramater that specifies the width of the DiracDelta function
potentialFunction = 'double-well';

%%%%%%%%%%%%%%%%% step 1: rough segmentation %%%%%%%%%%%%%%%%%%%%
%%% perform segmentation at each frame 
%%% by conventional method 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% step 2: contour initialization %%%%%%%%%%%%%%%%%%
%%% intialize contours in the first frame 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%% step 3: main loop (propagation) %%%%%%%%%%%%%%%%%
%%% propagate cell contours from frame K-1 to frame K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% step 3.1: conduct matching between frame {K, K+1, ..., K+Kt}


%%%%% step 3.2: build matching force map MF


%%%%% step 3.3: contours evolution 


%%%%% step 3.4: initialize contours for newly entered cells


%%%%% step 3.5: update information and prepare for next iteration 

