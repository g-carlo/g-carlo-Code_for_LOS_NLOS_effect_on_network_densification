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
% - single user positioned at the centre of the network (at the origin)
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

%%%%% PARAMETERS to be set

  % network geometry parameters

lambda      = 10/1e6;               % intestity of SPPP  (tested for values from 1 to 10^4 BSs/sqkm)
N_points    = 40e3;                 % number of simulated BSs
N_iteration = 10000;                % number of simulation snapshots  (suggested >= 10000)
  
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

    
%%%%%   Functions

PL_LOS = @(x)  K_L .* x .^( - beta_L);              % LOS path-loss function
PL_NLOS = @(x)  K_NL .* x .^( - beta_NL);           % NLOS path-loss function

p_LOS = @(x) exp( -(x/L).^2);                       % LOS likelihood function ~exp(-x^2)
p_LOS_exp = @(x) exp( -x/L);                        % LOS likelihood function ~exp(-x)
p_LOS_3GPP = @(x) 0.5-min(0.5.*ones(size(x)),5*exp(-156./x))+min(0.5.*ones(size(x)), 5*exp(-x./30));    % 3GPP LOS likelihood function


%%%%% Define variable to save data for all snapshots

SIR_vector = zeros(1,N_iteration);
rate_vector = zeros(1,N_iteration);
ASE_vector = zeros(1,N_iteration);


for n_iter = 1:N_iteration

        % generate network points
    numOfPoints = poissrnd(lambda.*area);                                                       % generate number  of points (i.e., BSs)
    position = edgeLength*( ( rand(1,numOfPoints) -0.5) + 1i*(rand(1,numOfPoints) - 0.5) );     % position of the points
    distance = abs(position);                                                                   % distance of the points from the origin
        
        % determine LOS likelihood depending on the distance
    switch LOS_function
        case 1
            p_of_having_LOS = p_LOS(distance);        % prob to have LOS                                                          
        case 2
            p_of_having_LOS = p_LOS_exp(distance);    % prob to have LOS                                                          
        case 3
            p_of_having_LOS = p_LOS_3GPP(distance);   % prob to have LOS                                                           
        otherwise 
            error('No other LOS likelihood functions');
    end


    drawVariable = rand(size(position));                % draw random Bernoulli variable for LOS/NLOS 
    LOS_boolean = drawVariable <= p_of_having_LOS;
    LOS_idx = find(LOS_boolean);                        % indices of LOS BSs
    NLOS_idx = find(~LOS_boolean);                      % indices of NLOS BSs
    P_rx = zeros(1,numOfPoints);                        % vector of user's received power from different BSs
    P_rx(LOS_idx)  = PL_LOS(distance(LOS_idx));         % user's received power from LOS BSs
    P_rx( NLOS_idx ) = PL_NLOS(distance(NLOS_idx));     % user's received power from NLOS BSs

    if numOfPoints > 0      
        [max_Prx, max_Prx_idx]= max(P_rx);              % compute the index of the serving BS (highest received power)

        interference_total = P_rx.* exprnd(fading_mean,1,length(P_rx));     % vector of total interference with added exponential fading component
        interference_total(max_Prx_idx) = 0;                                % remove serving BS from interferers
        interference_NLOS = interference_total(NLOS_idx);                   % LOS interference elements
        interference_LOS = interference_total(LOS_idx);                     % NLOS interference elements
        SIR_linear = max_Prx .* exprnd(fading_mean,1,1) ...                 % SIR of the user
            ./ sum( interference_total );    

        SIR_vector(n_iter) = SIR_linear;                                    % copy SIR for this snapshot into the SIR vector 
        rate_vector(n_iter) = log2(1+SIR_linear);                           % copy spectral efficiency (SE) for this snapshot into the SE vector 
        ASE_vector(n_iter) = numOfPoints*log2(1+SIR_linear)/area;           % copy Area Spectral Efficiency (ASE) for this snapshot into the ASE vector 
        
    end
    
    if fix(n_iter/N_iteration*100/10) > fix((n_iter-1)/N_iteration*100/10) 
        
        disp([ num2str(fix(n_iter/N_iteration*100/10)*10) '% completed'])  % display completed simulation perc.
        
    end

end

toc

disp(' ');
disp(['The coverage is     : ' num2str(100*sum(SIR_vector>SIR_threshold)/N_iteration) '%'] );
disp(' ');
disp(['The average rate is : ' num2str(mean(rate_vector)) ' bps'] );
disp(' ');
disp(['The ASE is          : ' num2str(mean( ASE_vector)*10^6) '*10^-6 bps/(Hz*sqm)'] );
disp(' ');


