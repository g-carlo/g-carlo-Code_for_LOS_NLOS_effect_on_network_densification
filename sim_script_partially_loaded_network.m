% This matlab script contains the code to test the stochatic-geometry based 
% model presented in the paper "Effect of LOS/NLOS Propagation on 5G Ultra-
% Dense Networks", submitted to "COMPUTER NETWORKS, Elsevier" and currently 
% under review.
% 
% The code provided in this script runs a simulation of the following
% SYSTEM MODEL:
% 
% - small-cell base stations deployed acccording to a homogeneous Spatial Poisson
% Point Process (SPPP) of density "lambda". The average number of BSs
% considered is "N_points"
% - we focus the analysis on a single user positioned at the centre of the
% network (at the origin)
% - we assume also an additional number of users, which are deployed
% aqccording to a homogeneous Spatial Poisson Point Process (SPPP) of
% density "lambda_users"
% - path-loss: dual-slope with LOS and NLOS paths, Rayleigh fading with
% exponentially distributed power ~exp(1)
% 
% The script returns:
% - Vector of SINR distribution of the user: "SIR_vector" 
% - Coverage as Prob[SIR > SIR_threshold]
% - Vector of Average spectral efficiency or rate of the typical user values: "rate_vector"
% - Vector of Area Spectral Efficiency (ASE) values: "ASE_vector"     

%%%%% Created by  :  Carlo Galiotto (galiotc@tcd.ie)
%%%%% Last update :  March 2017


%%%%% pre-initialization of script
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMTERS to be set  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % network geometry parameters

lambda       = 10/1e6;                % intestity of SPPP  (tested for values from 1 to 10^4 BSs/sqkm)
N_points     = 50000;                 % number of simulated BSs 
lambda_users = 1e3/1e6;               % density fo users (num of users is deterministic

N_iteration  = 10000;                 % number of simulation snapshots

  
  % channel model parameters

LOS_function    = 1;            % LOS likelihood function: 1 = ~exp(-x^2);    2 = ~exp(-x);    3 = 3GPP function
beta_L          = 2.09;         % LOS attenuation exponent   = 2.09
beta_NL         = 3.75;         % NLOS attenuation exponent  = 3.75
K_L_3GPP        = 103.8;        % LOS attenuation at 1 km
K_NL_3GPP       = 145.4;        % NLOS attenuation at 1 km
fading_mean     = 1;            % fading power (as long as fading is included, the fading mean does not impact the SINR)
L               = 69/sqrt(log(1/0.5));        % LOS likelihood parameter approx. 82.5m

  % PHY parameters

SIR_threshold   = 10^(-8/10);   % SIR threshold = -8dB

%%%%% Dependant PARAMETERS  

edgeLength = sqrt(  N_points / lambda);             % length of the square edge
area = edgeLength * edgeLength ;                    % area of the SPPP to simulate

K_L = 10^( -( K_L_3GPP - beta_L*10*3) / 10);        % LOS attenuation at 1m      
K_NL = 10^( -( K_NL_3GPP - beta_NL*10*3) / 10);     % NLOS attenuation at 1m 


tic                                                 % start timer

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PL_LOS = @(x)  K_L .* x .^( - beta_L);              % LOS path-loss function
PL_NLOS = @(x)  K_NL .* x .^( - beta_NL);           % NLOS path-loss function

p_LOS = @(x) exp( -(x/L).^2);                       % LOS likelihood function ~exp(-x^2)
p_LOS_exp = @(x) exp( -x/L);                        % LOS likelihood function ~exp(-x)
p_LOS_3GPP = @(x) 0.5-min(0.5.*ones(size(x)),5*exp(-156./x))+...
    min(0.5.*ones(size(x)), 5*exp(-x./30));         % 3GPP LOS likelihood function [36.814 - Urban pico-cells]


%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%
% Define variable to save data for all snapshots

SIR_vector = zeros(1,N_iteration);
rate_vector = zeros(1,N_iteration);
ASE_vector = zeros(1,N_iteration);


for n_iter = 1:N_iteration

    
    numOfBSs = poissrnd(lambda.*area);                      % generate number of BSs
    bs_position = edgeLength*( ( rand(1,numOfBSs) - 0.5) + 1i*(rand(1,numOfBSs) - 0.5) );           % position of the BSs
    distance_from_BS = abs( bs_position );                                                          % distance of the points from the origin
    
    % process only active BSs
    
    p_active =  (1 - (1 + ( lambda_users / (lambda * 3.5) ) )^-3.5  );      % probability of a BS being active
    drawVariable_active = rand(1,numOfBSs);                                 % draw random Bernoulli variable for active BSs
    active_BS_boolean = drawVariable_active < p_active;                     % find boolean indices of active BSs
    
            % determine LOS likelihood depending on the distance
    switch LOS_function
        case 1
            p_of_having_LOS = p_LOS(distance_from_BS);          % prob to have LOS                                                          
        case 2
            p_of_having_LOS = p_LOS_exp(distance_from_BS);      % prob to have LOS                                                          
        case 3
            p_of_having_LOS = p_LOS_3GPP(distance_from_BS);     % prob to have LOS                                                           
        otherwise 
            error('No other LOS likelihood functions');
    end

    drawVariable_LOS = rand(size(bs_position));                 % draw random Bernoulli variable for LOS/NLOS 
    LOS_boolean = drawVariable_LOS <= p_of_having_LOS;
    LOS_idx = find(LOS_boolean);                                % indices of LOS BSs
    NLOS_idx = find(~LOS_boolean);                              % indices of NLOS BSs
    P_rx = zeros(1,numOfBSs);                                   % vector of user's received power from different BSs
    P_rx(LOS_idx)  = PL_LOS(distance_from_BS(LOS_idx));         % user's received power from LOS BSs
    P_rx( NLOS_idx ) = PL_NLOS(distance_from_BS(NLOS_idx));     % user's received power from NLOS BSs
    
    if numOfBSs > 0      
        [max_Prx, serving_BSs_idx]= max(P_rx);                  % find index of serving base (the one providing the highest power)

        inactive_BSs_idx = find( ~active_BS_boolean );                          % find indices of inactive BSs
        interf_effective = P_rx.* exprnd(fading_mean,1,length(P_rx));           % vector of total interference with added exponential fading component
        interf_effective(unique([serving_BSs_idx inactive_BSs_idx]) ) = 0;      % remove serving BS and inactive BSs from interferers

        rx_pw_tmp = max_Prx .* exprnd(fading_mean,1,1);         % add exponential fading to received power from serving BS                            
        tot_int_tmp = sum( interf_effective );                  % total interference
        SIR_linear = rx_pw_tmp ./ tot_int_tmp;                  % compute SIR
                
        SIR_vector(n_iter) = SIR_linear;                                            % copy SIR for this snapshot into the SIR vector 
        rate_vector(n_iter) = log2(1+SIR_linear);                                   % copy spectral efficiency (SE) for this snapshot into the SE vector 
        ASE_vector(n_iter) = sum(active_BS_boolean)*log2(1+SIR_linear)/area;        % copy Area Spectral Efficiency (ASE) for this snapshot into the ASE vector         

    end
    
    if fix(n_iter/N_iteration*100/10) > fix((n_iter-1)/N_iteration*100/10) 
        
        disp([ num2str(fix(n_iter/N_iteration*100/10)*10) '% completed'])  % display completed simulation perc.
        
    end

end 

toc

% display results
disp(' ');
disp(['The coverage is     : ' num2str(100*sum(SIR_vector>SIR_threshold)/N_iteration) '%'] );
disp(' ');
disp(['The average rate is : ' num2str(mean(rate_vector)) ' bps'] );
disp(' ');
disp(['The ASE is          : ' num2str(mean( ASE_vector)*10^6) '*10^-6 bps/(Hz*sqm)'] );
disp(' ');













