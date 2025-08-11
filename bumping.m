%estimating delta using MC and central differences:
clear all;
format long;
K = 100;
T = 1;
vol = 0.2;
r = 0.05;
t0 = 0;
S0 = 98;
delta = (cblackscholes(K,T,0,r,vol,S0+0.01,t0)-cblackscholes(K,T,0,r,vol,S0-0.01,t0))/(2*0.01);
for k=1:8
    N(k) = 4^k;
    ds(k) = 2^(-k);
end
%Naive scheme: 
for i=1:8
    %for each choice of sample size N and stock step ds
    for j=1:8
        estimate_naive(i,j) = ApproxDelta_naive(K,T,vol,r,t0,S0,N(j),ds(i));
        error_naive(i,j) = estimate_naive(i,j)-delta; %compared with "exact" delta
    end
    h(i) = plot(log(N), error_naive(i,:),'--o'); %plot of errors (error vs ds for different n)s
    hold on;
end
title("Plot of errors using naive scheme");
xlabel("log(N)");
ylabel("error");
grid on;
% Create legend with step size labels
legendStrings = arrayfun(@(d) sprintf("ds = %.4f", d), ds(1:8), 'UniformOutput', false);
legend(h, legendStrings);

disp("Giles test_")
ApproxDelta_Giles(K,T,vol,r,t0,S0,100,0.1)

figure
%Giles scheme :
for i=1:8
    %for each choice of sample size N and stock step ds
    for j=1:8
        estimate_giles(i,j) = ApproxDelta_Giles(K,T,vol,r,t0,S0,N(j),ds(i));
        error_giles(i,j) = estimate_giles(i,j) -delta; %compared with "exact" delta
    end
    h(i) = plot(log(N), error_giles(i,:),'--o'); %plot of errors (error vs ds for different n)s
    hold on;
end
title("Plot of errors using Gile's scheme");
xlabel("log(N)");
ylabel("error");
grid on;
% Create legend with step size labels
legendStrings = arrayfun(@(d) sprintf("ds = %.4f", d), ds(1:8), 'UniformOutput', false);
legend(h, legendStrings);


%% NUMMERICAL CONCLUSIONS:
%% Although both methods are convergent, we see a clear flaw from the first one (naive) which is that the estimator used
%% is very inefficient in the sense that the variance of the estimator is of the order of 1/N*ds^2 which means that as ds gets small
%% the sample size needed (N) to achieve a prescribed level of convergence explodes making it hopeless to be convergent in practice probably.
%% However, this is as discussed by Giles due to the fact that the realisations used to average out the payoff at maturity are independed from one another,
%% if the same realisations of the brownian motion are used to create realisations of the up and down initial prices perturbations, then they will give a better estimator since the errors will "cancel out"
%% Therefore, while the naive seems simpler to implement and more intuitive at first, it turns out that as ds -> 0 (N fixed) the error does not necessarily go down
%% while for the Giles approach, the error goes down.













function estimated_delta = ApproxDelta_naive(K,T,vol,r,t0,S0,n,ds) 
    %computes the approx delta using MC for the value and center diff for
    %the derivative (2 error sources essentially, one random and the other
    %deterministic)
    estimated_delta = 0;
    right = PriceEuroCallMC(K,T,vol,r,t0,S0+ds,n);
    left = PriceEuroCallMC(K,T,vol,r,t0,S0-ds,n);
    estimated_delta = (right-left)/(2*ds);
end


function estimated_delta = ApproxDelta_Giles(K,T,vol,r,t0,S0,n,ds)
    %computes the approx delta by using the same stream of random numbers
    %for both simulations of the price
    common_realisations = randn(n,1);
    right_final= (S0+ds)*exp((r-vol^2/2)*(T-t0) + vol*sqrt(T-t0)*common_realisations); %for the right
    score_right= max(right_final-K,zeros(n,1)); %payoff
    V_right= exp(-r*(T-t0))*mean(score_right); %discounted ,average payoff
    left_final = (S0-ds)*exp((r-vol^2/2)*(T-t0) + vol*sqrt(T-t0)*common_realisations); %for the left 
    score_left= max(left_final-K,zeros(n,1)); %payoff
    V_left = exp(-r*(T-t0))*mean(score_left); %discounted ,average payoff
    estimated_delta = (V_right-V_left)/(2*ds);
    
end
