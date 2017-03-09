% This matlab function contains the code to test the stochatic-geometry based 
% model presented in https://arxiv.org/abs/1507.01757
% 
% The code provided in this script runs a simulation of the following
% SYSTEM MODEL:
% 
% - small-cell base stations deployed acccording to a homogeneous Spatial Poisson
%   Point Process (SPPP) of density "lambda_B". The average number of BSs
%   considered is "N_points"
% - we assume users are deployed according to a homogeneous Spatial Poisson 
%   Point Process (SPPP) of density "lambda_U"
% - users connect to the BSs from which the received power is the strongest
%
% The function computes the probability of a base station being active
% "p_active", i.e., the probability of a base station to have at least one
% user to serve. The computation is based on the model descrived above.
%  
% The function inputs are:
%   1) "lambda_B"   : density of base stations
%   2) "lambda_U"   : density of users
%   3) "N_points"   : number of base stations to simulate,
%   4) "n_snapshot" : number of snapshot per couple (lambda_B,lambda_U)
%
% The function outputs are:
%   1) "p_active_mean"  : matrix ( number of BS density values x number of BS density values)
%                         with average probability of a base station being
%                         active "p_active"; one value per couple (lambda_B,lambda_U)
%   2) "p_active_std"   : matrix ( number of BS density values x number of BS density values)
%                         with standard deviation of the probability of a base station being
%                         active "p_active"; one value per couple (lambda_B,lambda_U)

%%%%% Created by  :  Carlo Galiotto (galiotc@tcd.ie)
%%%%% Last update :  March 2017


function [p_active_mean, p_active_std] = fnc_check_p_Active_LOS_NLOS(lambda_B,lambda_U,N_points,n_snapshot)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMTERS to be set  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is meant for the following values

perc_of_interest        = 0.8/2;           % percentage (80%) of the perimeter used for collecting results, 20% discared to avoid edge effect
    
beta_L = 2.09;                      % LOS attenuation exponent   = 2.09
beta_NL = 3.75;                     % NLOS attenuation exponent  = 3.75    

K_L_3GPP = 103.8;                   % LOS attenuation at 1 km
K_NL_3GPP = 145.4;                  % NLOS attenuation at 1 km


%%% DEPENDANT PARAMETERS  

edgeLength = sqrt(  N_points / lambda_B);  % length of the square edge
area = edgeLength * edgeLength ;           % area of the SPPP to simulate

    % channel model

K_L = 10^( -( K_L_3GPP - beta_L*10*3) / 10);            % LOS attenuation at 1m      
K_NL = 10^( -( K_NL_3GPP - beta_NL*10*3) / 10);         % LOS attenuation at 1m      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PL_LOS = @(x)  K_L .* x .^( - beta_L);                  % LOS path-loss function
PL_NLOS = @(x)  K_NL .* x .^( - beta_NL);               % NLOS path-loss function

p_LOS_3GPP = @(x) 0.5-min(0.5.*ones(size(x)),5*exp(-156./x))+...
    min(0.5.*ones(size(x)), 5*exp(-x./30));             % 3GPP LOS likelihood function [36.814 - Urban pico-cells]


%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%
% Define variable to save data for all snapshots

p_active_sim_vect = zeros( 1, n_snapshot);

for snapshot_idx = 1:n_snapshot
    
    %% deploy base stations and users
    numOfBSs = poissrnd(lambda_B.*area);            % generate number of BSs
    numOfUEs = poissrnd(lambda_U.*area);            % generate number of UEs
    
    BS_of_interest = false(1,numOfBSs);             % vector of the BSs to be included in the results
    BS_to_UE_mapping = false(numOfBSs,numOfUEs);    % create (empty) BS to UE association matrix (which UE is connected to which BS)
%     UE_of_interest = false(1,numOfUEs);

    position_BS = edgeLength*( ( rand(1,numOfBSs) -0.5) + 1i*(rand(1,numOfBSs) - 0.5) );         % position of the BSs
    position_UE = edgeLength*( ( rand(1,numOfUEs) -0.5) + 1i*(rand(1,numOfUEs) - 0.5) );         % position of the UEs
    
    %% find base stations within the area of interest 
    % for each BS, check which BSs are contained within the area of interest
    % (some BSs will be ignored to remove edge-effects)
    for k = 1:numOfBSs
        BS_of_interest(k) = inpolygon( real( position_BS(k) ) , imag( position_BS(k) ), ...
            [ -perc_of_interest*edgeLength perc_of_interest*edgeLength] , [ -perc_of_interest*edgeLength perc_of_interest*edgeLength] ); 
    end
    
    %% compute distance from BS to UE
    tmp = position_BS.';    
    pos_BS_mtx_MxK = tmp( :, ones( 1,numOfUEs ) );          % create a matrix numOfBSs x numOfUsers with positions (repetitions) of BSs
    tmp = position_UE;
    pos_UE_mtx_MxK = tmp( ones(1,numOfBSs), :);             % create a matrix numOfBSs x numOfUsers with positions (repetitions) of UEs
    
    distance_BS_to_UE = abs( pos_BS_mtx_MxK - pos_UE_mtx_MxK );     % compute the distace from each BS to each UE

    %% compute signal power from BS to UE
    p_of_having_LOS = p_LOS_3GPP(distance_BS_to_UE);                      % determine LOS likelihood depending on the distance
    drawVariable = rand(size(distance_BS_to_UE));                         % draw random Bernoulli variable for LOS/NLOS   
    LOS_boolean = drawVariable <= p_of_having_LOS;                        
    P_rx = zeros(numOfBSs,numOfUEs);                                      % create empty matrix of user's received power from each BS
    P_rx( LOS_boolean  )  = PL_LOS(distance_BS_to_UE( LOS_boolean ));     % received power from BSs in LOS 
    P_rx( ~LOS_boolean ) = PL_NLOS(distance_BS_to_UE( ~LOS_boolean ));    % received power from BSs in NLOS

    %% compute the element of the BS to UE association matrix 
    [ ~, tmp_max_pw_idx ] = max( P_rx );                                  % find index of serving base (the one providing the highest power)
    BS_to_UE_mapping( sub2ind( [numOfBSs, numOfUEs] , tmp_max_pw_idx, 1:numOfUEs ) ) = true;        % write "1" on empty matrix to associate UE and serving BS
    BS_with_UEs = sum( BS_to_UE_mapping,2 );                              % find number of users served by each BS
    BS_active = BS_with_UEs > 0 ;                                         % find the active BSs, i.e., the BSs serving at least 1 UE
    
    %% compute the probability of a base station being active "p_active"
    % compute "p_active" and save value for this snapshot
    p_active_sim_vect( snapshot_idx )  = sum( BS_active( BS_of_interest ) ) / sum( BS_of_interest  );   
   
    
end

% compute average "p_active" (averaged over snapshots)
p_active_mean = mean( p_active_sim_vect);

% compute standard deviation of "p_active" (averaged over snapshots)
p_active_std = std( p_active_sim_vect );





