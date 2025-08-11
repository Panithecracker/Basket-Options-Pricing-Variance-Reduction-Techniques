%% Studying weak convergence rate of Euler Maruyama method to price Vanilla calls
%% in contrast to barrier options on an OU underlying. Also implementing Maruyama combined with 
%% Brownian bridge and barrier shiftings to see how the converge rate is increased from 
%% order 0.5 to 1 too.

clear all 
format long

%% Underlying parameters:
%mean reversion time, long term mean, volatility, barrier, current price
%, strike, time till expiration and risk free rate:
S0 = 14;
T = 2;
K = 14;
B0 = 13.5;
sigma = 0.5;
kappa = 2;
theta = 14;
r = 0;
M = [8,16,32,64,128,256]; %partition orders



% %% VANILLA OPTION OF OU 
% for i=1:length(M)
%     [V,ster,CPUt,varsc,eb,h] = ouscratch(1e6,M(i),-1,1,false);
%     H(i,1) = log2(h);
%     bias_error(i,1) = log2(eb); 
%     Value(i,1) = log2(V(1,1));
% end
% 
% plot(H,bias_error);
% hold on
% plot([H(5,1) H(1,1)] ,[bias_error(5,1) bias_error(5,1)+1*4],'--');
% axis equal
% grid on
% xlabel("log2(h)");
% ylabel("log2(bias)");
% title("Bias error estimate vanilla call option");
% 
% figure
% plot(H,Value);
% axis equal;
% grid on;
% title("Price estimate vanilla call option");
% xlabel("log2(h)");
% ylabel("log2(V(h))");





%% DOWN-OUT BARRIER OPTION OF OU




%% Now, we use Euler Maruyama + Bbdrige technique:
%For this case we have seen in the theory that this combination leads to a
%1 order method better than the 0.5 before.
%We will simulate OU at T (known mean and variance) and use the sample mean
%as the exact solution to compare with the Bbride + Maruyama accuracy and
%convergence rate
N = 1000000;
for i=1:length(M) %different time step
    dt= T/M(i); 
    H(i,1) = log2(dt);
%     S= NaN(N,M(i)+1); 
    S = zeros(N,2);
    S(:,1)= S0;
    B= B0; 
%     B= B0-0.5826*sigma*sqrt(dt); %upwards in this problem 
    %% EULER MARUYAMA WITH BROWNIAN BRIDGE 
    product = ones(N,1); %compound product from brownian bridge method discount factor
%     dW = randn(N,M(i)); %wiener increments matrix (each row is a path)
    for j=1:M(i)
        %% maruyama+bbridge using only N*2 matrix (current and next)
        dW = randn(N,1);
        S(:,2) = S(:,1) + kappa*(theta-S(:,1))*dt + sigma*sqrt(dt)*dW(:,1);
        product = product .*(1-exp((-2*max(0,S(:,2)-B).* max(0,S(:,1)-B))/(sigma*sigma*dt)));
        S(:,1) = S(:,2);
%         S(:,j+1)= S(:,j) + kappa*(theta-S(:,j))*dt + sigma*sqrt(dt)*dW(:,1);
%         t = j*dt;
%         std =  sqrt((sigma^2*(1-exp(-2*kappa*t)))/(2*kappa));
%         S(:,j+1) = (S0*exp(-kappa*t)+theta*(1-exp(-kappa*t)))*ones(N,1)+ std*randn(N,1);
%         product = product .* (1-exp((-2*max(0,S(:,j+1)-B0).* max(0,S(:,j)-B0))/(sigma*sigma*dt)));
    end
    payoff= zeros(N,1);
    payoff= exp(-r*T)*max(0,S(:,end)-K).*product; %discounted payoffs 
    V(i)= mean(payoff); %option price estimate
    varsc(i)= var(payoff); 
    ster(i)= 3*sqrt(varsc(i)/N); %standard error (99 % sure not more than this)
    %% estimate bias error:
    if i > 1
        bias_error(i-1,1) = log2(2*(V(i)-V(i-1)));
    end
    
end

plot(H(1:5,1),bias_error);
hold on
plot([H(5,1) H(1,1)] ,[bias_error(5,1) bias_error(5,1)+1*4],'--');
axis equal
grid on
xlabel("log2(h)");
ylabel("log2(bias)");
title("Bias error estimate barrier option using Brownian Bridge");



% 
% conv_rate = zeros(length(M),1);
% conv_prec = zeros(length(M),1);
% for j=3:length(M)
%    conv_rate(j,1) = log2((V(j-1)-V(j-2))/(V(j)-V(j-1)));
%    conv_prec(j,1) = log2((Preciomedio(j-1)-Preciomedio(j-2))/(Preciomedio(j)-Preciomedio(j-1)));
% end
% C(programs,1) = conv_rate(3,1);
% Cprec(programs,1) = conv_prec(3,1);

% hist(C);

% figure
% plot(H,bias_error);
% hold on
% plot([H(5,1) H(1,1)] ,[bias_error(5,1) bias_error(5,1)+0.5*4],'--');
% axis equal
% grid on
% xlabel("log2(h)");
% ylabel("log2(bias)");
% title("Bias error estimate barrier call option");
% 
% figure
% plot(H,Value);
% axis equal;
% grid on;
% title("Price estimates barrier call option");
% xlabel("log2(h)");
% ylabel("log2(V(h))");
% 
