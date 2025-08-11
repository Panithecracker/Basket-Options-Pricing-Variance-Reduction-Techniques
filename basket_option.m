%monte carlo simulation to price a basket option consisting of a
%particulatr fraction of n stocks (specified by weights) at a strike of K
%with a particular initial prices and based on geometric brownian motion
%that includes stocks volatilities, correlation matrix and the risk free rate r

format long 
clear all

%% the exact solution is the expected payoff of the basket discounted by r

%We study two different basket types:

%%  BASKET (n) parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=1; %time to expiration
num_stocks = 10; 
r = 0.1;
K = 100;
rho = 0.5; %correlation
corr = rho * ones(num_stocks,num_stocks) + (1-rho) * eye(num_stocks); %correlation matrix
weights = (1/num_stocks)*ones(1,num_stocks); %weight for the underlyings
S0 = 100*ones(num_stocks,1); %starting prices of the underlyings
sigma = 0.2*ones(num_stocks,1); %volatilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TECH BASKET PARAMETERS:
% T = 1;
% num_stocks = 4; 
% r=0.01; %risk free rate
% K=23; %basket's strike price
% corr = [1 0.55 0.53 0.51; 0.55 1 0.55 0.48; 0.53 0.55 1 0.47; 0.51 0.48 0.47 1];
% weights = (1/num_stocks)*ones(1,num_stocks); %weight for the underlyings
% S0 = [25.87;26.77;24.54;18.63];
% sigma = [0.204;0.207;0.211;0.258];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATED GBM PARAMETERS: 
n=100; %partition order
dt= T/n; %time step
V = chol(diag(sigma)*corr*diag(sigma),'lower'); %cholesky decomposition to simulate correlated brownians
S(:,1) = S0;



% %extra: animation
% figure;
% hold on;
% axis([0,1,-10,30]);
% for s=1:num_stocks
%     p_h(s,1) = scatter(0,0,'filled'); %offset from start price 
% end


%% CONTROL VARIATE SETUP : FINDING OPTIMAL LAMBDAS
%solving approximate linear system by approximating covariance matrix for control variate:
mean_maturity = S0*exp(r*T); %mean stock prices at maturity
% tsetup = tic;  % Start timer for MMC set up 
 for i=1:5000
    for j=1:n
        S = S + r*S*dt + sqrt(dt)*S.*V*randn(num_stocks,1);
        S = max(S,0);
    end
    Payoffs(i,1) = max(weights*S-K,0);

    for p=1:num_stocks
        Temp = mean_maturity; 
        Temp(p,1) = S(p,1); 
        Perizzellis(i,p) = max(weights*Temp-K,0);
    end
    S(:,1) = S0;
 end
 %estimate covariances
 M = cov(Perizzellis);
 b = zeros(num_stocks,1);
for p=1:num_stocks
    cov_value = (-1)*cov(Payoffs, Perizzellis(:,p));
    b(p,1) = cov_value(1,2);
end
lambda = linsolve(M,b); %solution of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1: SIMULATE SAMPLE MEANS OF THE PAYOFF (CRUDE MC) AND OF THE MEAN MC FOR 100 TIMES WITH DIFFERENT N:        
%% MEASURE TIME OF EXECUTION FOR EACH METHOD AND ESTIMATE THE CORRELATION BETWEEN PAYOFF AND THE TOTAL CONTROL VARIATE 

sample_sizes = 1000:1000:10000; 
results = zeros(length(sample_sizes), 5); 
mean_maturity = S0*exp(r*T); %mean stock prices at maturity
%calculate the mean of each Perizzelli control variate using Black Scholes
%formula for a vanilla call option
cv_means = zeros(num_stocks,1);
 for p=1:num_stocks
  ficticious_k = (K + weights(1,p)*mean_maturity(p,1) - weights*mean_maturity) / (weights(1,p));
  cv_means(p,1) = weights(1,p)*exp(r*T)*blsprice(S0(p,1), ficticious_k, r, T, sigma(p,1));
 end
%  time_setup = toc(tsetup); %time taken to set up the control variate estimator
%  MMC_Time = time_setup;
%  MC_Time = 0;

%% MAIN PRICING SIMULATION:
for s = 1:length(sample_sizes)
    N = sample_sizes(s); 
    
    Payoffs = zeros(N,1); %observation payoffs
    MMC = zeros(N,1);  %observation of the estimator induced from combination of control variates

    for i=1:N %simulate N basket outcomes
%         sim_timer = tic;  % Start timing the simulation
        for j=1:n
            S = S + r*S*dt + sqrt(dt)*S.*V*randn(num_stocks,1);
            S = max(S,0);
        end
        Payoffs(i,1) = max(weights*S-K,0); %compute payoff of this outcome
%         MC_Time = MC_Time + toc(sim_timer);
        %compute observed control variates
        MMC(i,1) = Payoffs(i,1);
        for p=1:num_stocks
            Temp = mean_maturity; 
            Temp(p,1) = S(p,1); 
            MMC(i,1) = MMC(i,1)+ lambda(p,1)*(max(weights*Temp-K,0)-cv_means(p,1));
        end
        S(:,1) = S0; %reset initial prices
%         MMC_Time = MMC_Time + toc(sim_timer);  % Include lambda time
    end
    %compute the estimate for the mean and variance using both MC and MMC
    %estimators:
    MC_mean = exp(-r*T)*mean(Payoffs);
    MC_var = (exp(-2*r*T)*var(Payoffs))/N;
    MMC_mean = exp(-r*T)*mean(MMC);
    MMC_var = (exp(-2*r*T)*var(MMC))/N;
    var_reduction = ceil(100*(MMC_var-MC_var)/MC_var); %percent of variance reduction

    % Store results for sample size N in table
    results(s, :) = [MC_mean, sqrt(MC_var), MMC_mean,sqrt(MMC_var),var_reduction];
%     execution_times(s,:) = [MC_Time, MMC_Time];
%     MC_Time = 0;
%     MMC_Time = time_setup;
    [a, p] = corrcoef(Payoffs, MMC-Payoffs);
    kappa(s) = a(1, 2);
    fprintf("Done with N = %d\n",N);
end

% Convert to table for better readability
results_table = array2table(results, 'VariableNames', {'MC', 'MC_deviation', 'MMC', 'MMC_deviation','Var reduction (%)'}, 'RowNames', string(sample_sizes));
disp(results_table);

%plot of Cpu times for each method:
% figure
% title("Computational costs relative to the sample size");
% plot([1:10],execution_times(:,1)');
% hold on;
% plot([1:10],execution_times(:,2)'-time_setup);
% legend("MC","MMC");
% xlabel("N/10");
% ylabel("TIME COST");
% grid on;
% 



%% 3 STUDYING THE VARIANCE REDUCTiON WITH VARYING BASKET SIZE (n):
%computing several sample means of size N= 5000


N = 1000;
S(:,1) = S0; %reset initial prices
Payoffs = zeros(N,1); %observation payoffs
MMC = zeros(N,1);  %observation of the estimator induced from combination of control variates
for t = 1:500
    for i=1:N %simulate N basket outcomes
        for j=1:n
            S = S + r*S*dt + sqrt(dt)*S.*V*randn(num_stocks,1);
            S = max(S,0);
        end
        Payoffs(i,1) = max(weights*S-K,0); %compute payoff of this outcome
        %compute observed control variates
        MMC(i,1) = Payoffs(i,1);
        for p=1:num_stocks
            Temp = mean_maturity; 
            Temp(p,1) = S(p,1); 
            MMC(i,1) = MMC(i,1)+ lambda(p,1)*(max(weights*Temp-K,0)-cv_means(p,1));
        end
        S(:,1) = S0; %reset initial prices
    end
    %compute the estimate for the mean and variance using both MC and MMC
    %estimators:
    MC_means(t,1) = exp(-r*T)*mean(Payoffs);
    MC_vars(t,1) = (exp(-2*r*T)*var(Payoffs))/N;
    MMC_means(t,1) = exp(-r*T)*mean(MMC);
    MMC_vars(t,1) = (exp(-2*r*T)*var(MMC))/N;
    


end

% Compute histogram data
figure;
hold on;

% Plot histograms (with transparency for clarity)
histogram(MC_means, 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'EdgeColor', 'black');
histogram(MMC_means, 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'black');
% Fit normal distributions
pd_MC = fitdist(MC_means, 'Normal'); % Fit normal to raw Monte Carlo samples
pd_MMC = fitdist(MMC_means, 'Normal'); % Fit normal to Control Variate samples
% Generate smooth x values
x_values = linspace(min([MC_means; MMC_means]), max([MC_means; MMC_means]), 100);
% Compute normal PDF for both distributions
y_MC = normpdf(x_values, pd_MC.mu, pd_MC.sigma);
y_MMC = normpdf(x_values, pd_MMC.mu, pd_MMC.sigma);
% Overlay fitted normal curves
plot(x_values, y_MC, 'b', 'LineWidth', 2);   % Blue curve for MC_means
plot(x_values, y_MMC, 'r', 'LineWidth', 2);   % Red curve for MMC_means
% Labels and title
xlabel('Value');
ylabel('Probability Density');
title('Variance Reduction illustration for N = 1000');
legend({'MC estimator', 'Pezirrelli estimator'}, 'Location', 'best');
hold off;













