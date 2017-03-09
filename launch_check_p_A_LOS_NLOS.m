% This matlab script contains the code to test the stochatic-geometry based 
% model presented in the paper "Effect of LOS/NLOS Propagation on 5G Ultra-
% Dense Networks", submitted to "COMPUTER NETWORKS, Elsevier" and currently 
% under review.
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
% The script aims at:
% - Comparing the analytical expression of the probability of a base station being active, i.e., 
%   p_active =  (1 - (1 + ( lambda_U / (lambda_B * 3.5) ) )^-3.5  ) with the simulated results for 
%   the same probability with LOS/NLOS propagation
% - Plotting the analytical expression of "p_active" against the ratio lambda_U/lambda_B;
% - Plotting the simulated curves of of "p_active" for different combinations 
%   of (lambda_B,lambda_U) against the ratio lambda_U/lambda_B;


%%%%% Created by  :  Carlo Galiotto (galiotc@tcd.ie)
%%%%% Last update :  March 2017


%%%%% pre-initialization of script
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMTERS to be set  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_snapshot              = 100;                   % number of snapshots  (NOTE, it takes long!)

lambda_B_vector = 10.^( 1:0.2:4 ) /1e6;          % vector of BS density values in BSs per sqkm
lambda_U_vector_over_lambda_B = 10.^( -1.6:0.2:1 );         % vector of (relative) UE density values: this value will be 
                                                            %    multiplied by to BS density to obtain the UE density, i.e.,
                                                            %    "actual UE density" =  "lambda_B" x "lambda_U" 
N_points = [ 7e3 7e3 7e3 7e3 5e3 ...             % number of BSs to simulate (this number depends 
             5e3 5e3 5e3 5e3 3e3 ...             %    on the density of BSs and UEs, and it is such that
             3e3 2e3 1e3 1e3];                   %    to limit the computational complexity.
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

p_active_fnc = @(lam_U,lam_B)  (1 - (1 + ( lam_U ./ (lam_B * 3.5) ) ).^-3.5  );   % probability of a BS being active
 

%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE VARIABLES   %%%%%%%%%%%%%%%%%%%%%%%%%%
% Define variables where to store the values of probabilities p_A

p_A_mean_matrix = zeros( length(lambda_B_vector ),length( lambda_U_vector_over_lambda_B) );
p_A_std_matrix = zeros( length(lambda_B_vector ),length( lambda_U_vector_over_lambda_B) );

cnt = 1;        % set counter to 1 (to diplay simulation progress)

tic             % start timer 
for lam_B_idx = 1:length(lambda_B_vector)           % iterate over BS densities
    
    for lam_U_idx = 1:length(lambda_U_vector_over_lambda_B)       % iterate over UE densities
    
        lambda_B = lambda_B_vector( lam_B_idx );    
        lambda_U = lambda_U_vector_over_lambda_B( lam_U_idx );
        n_point_current = N_points( lam_U_idx );
        % call function that computes p_active
        [p_active_mean, p_active_std] = fnc_check_p_Active_LOS_NLOS(lambda_B,lambda_B*lambda_U,n_point_current,n_snapshot);
        p_A_mean_matrix( lam_B_idx, lam_U_idx ) = p_active_mean;    % save mean value of p_A_mean_matrix (for given BS and UE densities)
        p_A_std_matrix( lam_B_idx, lam_U_idx ) = p_active_std;      % save standard deviation of p_A_mean_matrix (for given BS and UE densities)
        
        cnt = cnt + 1;                                              % increment counter (to diplay simulation progress)
        
        % display completed simulation perc.
        N_iteration = length(lambda_B_vector ) * length( lambda_U_vector_over_lambda_B);
        if fix(cnt/N_iteration*100/10) > fix((cnt-1)/N_iteration*100/10) 
        	disp([ num2str(fix(cnt/N_iteration*100/10)*10) '% completed'])  
        end
        
    end
    
end
toc             % stop timer 


%%%%%%%%%%%%%%%%%%%%%%%%  PLOT RESULTS   %%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
fig_handles.h(1) = 0;

for n_row =  1:size(p_A_mean_matrix,1)
  
    % for each value of lambda_B, plot p_active against lambda_U/lambda_B
    fig_handles.h(n_row) = semilogx(lambda_U_vector_over_lambda_B, p_A_mean_matrix(n_row,:),'linewidth',1);
    hold on;
    
end

% plot the analytical expression of p_active against lambda_U/lambda_B
fig_handles.h(n_row+1) = semilogx(lambda_U_vector_over_lambda_B, p_active_fnc( lambda_U_vector_over_lambda_B,1 ) ,'-g+','linewidth',1);
grid on
xlabel('\lambda_U/\lambda_B');
ylabel('Probability p_A');
legend([fig_handles.h(1) fig_handles.h(n_row+1) ],{'Simulation - LOS/NLOS model','Analytical model - single slope'},'location','SouthEast');



% for n_row =  1:size(p_A_mean_matrix,1)
%    
%     semilogx(lambda_U_vector_over_lambda_B, abs( p_A_mean_matrix(n_row,:) - p_active_fnc( lambda_U_vector_over_lambda_B,1 ) ) )
%     hold on;
%     
% end

% generate a small-sub-plot to zoom into the curve

fig1 = figure(1);
fig_small = axes('position',[ .18 .51 .4 .4]);      % create new pair of axes
set(fig_small,'FontSize',8);
box on                                              % put box around new pair of axes
indexOfInterest = 10:12;                            % indices of lambda_U of interest
for n_row =  1:size(p_A_mean_matrix)
   
    % for each value of lambda_B, plot p_active against lambda_U/lambda_B
    % (only for indices of lambda_U indicated above)
    semilogx(lambda_U_vector_over_lambda_B(indexOfInterest), p_A_mean_matrix(n_row,indexOfInterest),'linewidth',1)
    hold on;
    
end

% plot the analytical expression of p_active against lambda_U/lambda_B
semilogx(lambda_U_vector_over_lambda_B(indexOfInterest), p_active_fnc( lambda_U_vector_over_lambda_B(indexOfInterest),1 ) ,'-g+','linewidth',1)
axis tight


% Create textbox
annotation(fig1,'textbox',...
    [0.333142857142857 0.604761904761906 0.091857142857143 0.0714285714285756],...
    'String',{'0.0192'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create arrow
annotation(fig1,'arrow',[0.380357142857142 0.380357142857142],...
    [0.66666666666667 0.726190476190486],'HeadStyle','plain');

% Create arrow
annotation(fig1,'arrow',[0.380357142857142 0.380357142857142],...
    [0.825238095238111 0.765714285714305],'HeadStyle','plain');

% Create line
annotation(fig1,'line',[0.380357142857142 0.380357142857142],...
    [0.837095238095241 0.664285714285717]);
