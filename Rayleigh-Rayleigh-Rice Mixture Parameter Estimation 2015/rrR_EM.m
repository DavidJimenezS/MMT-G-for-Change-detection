%This code was made from the article 
% "A theoretical framework for change detection based on a compound 
% multiclass statistical model of the difference image "

% if you used remember to cited bibtex here:
% @article{zanetti2017theoretical,
%   title={A theoretical framework for change detection based on a compound multiclass statistical model of the difference image},
%   author={Zanetti, Massimo and Bruzzone, Lorenzo},
%   journal={IEEE Transactions on Geoscience and Remote Sensing},
%   volume={56},
%   number={2},
%   pages={1129--1143},
%   year={2017},
%   publisher={IEEE}
% }
%%
clear all
close all
clc

%% Initialization

%________________ lake data ___________
load('mulgaria_95_96.mat')
%load('fire_2013.mat')

%--------lake & fire ------------
I = lake';
Aspot = lakef';
clear lake lakef;

%________________ contest data ___________
%load('data_contest.mat')
%load('mulargia_2001_2003.mat')

%_________________Magnitude of difference _______________


ro = sqrt((I(:) - Aspot(:)).^2);

ro(ro == 0) = eps;


N = length(ro);



%% Population of W11, W12 and W2 for ML initial parameters

%Qunatile thresholding
nbins = 300;
[counts ,centers] = hist(ro,nbins);
counts_n = counts./N;
Y = quantile(ro,counts_n);
y_d = diff(Y)./diff(counts_n);

T = centers(y_d == max(y_d));

idx_w1 = ro <= T;
idx_w2 = ro > T;
W1 = ro(idx_w1);
W2 = ro(idx_w2);

change_map = ro;
change_map(idx_w1) = 0;
change_map(idx_w2) = 1;


%median thresholding for multi-class no change
T = (max(W1) - min(W1))/2;

idx_w11 = W1 <= T;
idx_w12 = W1 > T;
W11 = W1(idx_w11);
W12 = W1(idx_w12);

[m_c n_c] = size(Aspot);

figure, imshow(reshape(ro,m_c,n_c)), colorbar
figure, imshow(reshape(change_map,m_c,n_c)), colorbar

alpha_1 = length(W11)/N;
alpha_2 = length(W12)/N;
delta_1 = raylfit(W11);
delta_2 = raylfit(W12);

pd = fitdist(W2,'Rician');
nu = pd.s;
delta = pd.sigma;

%% Compute log-likelihood

p_x_w11 = rayleigh(ro,delta_1);
p_x_w12 = rayleigh(ro,delta_2);
p_x_w2 = rician(ro,delta,nu);

%Solve problem with low probabilities round to zero
p_x_w11(p_x_w11 == 0) = eps;
p_x_w12(p_x_w12 == 0) = eps;
p_x_w2(p_x_w2 == 0) = eps;
p_x_w2(isnan(p_x_w2)) = eps;
p_x_w2(isinf(p_x_w2)) = eps;


p_ro_Psi = (alpha_1*p_x_w11) + ...
    (alpha_2*p_x_w12) + ...
    ((1 - alpha_1 - alpha_2)*p_x_w2);

clear idx_w1 idx_w2

%imagesc(reshape(p_ro_Psi,m_c,n_c))


log_likelihood(1) = sum(log(p_ro_Psi));

log_likeliold = 100;

figure;
hold on

Q_psi_psihat = zeros(N,3);
k = 1;

while abs(log_likelihood(k) - log_likeliold)/log_likeliold > 10e-8
    
    
    
    %%Expectation step
    p_x_w11 = rayleigh(ro,delta_1);
    p_x_w12 = rayleigh(ro,delta_2);
    p_x_w2 = rician(ro,delta,nu);
    
    %Solve problem with low probabilities round to zero
    p_x_w11(p_x_w11 == 0) = eps;
    p_x_w12(p_x_w12 == 0) = eps;
    p_x_w2(p_x_w2 == 0) = eps;
    p_x_w2(isnan(p_x_w2)) = eps;
    p_x_w2(isinf(p_x_w2)) = eps;
    
    
    p_w11_x_psi = (alpha_1*p_x_w11)./p_ro_Psi;
    p_w12_x_psi = (alpha_2*p_x_w12)./p_ro_Psi;
    p_w2_x_psi = ((1 - alpha_1 - alpha_2)*p_x_w2)./p_ro_Psi;
    
    Q_psi_psihat(:,1) = p_w11_x_psi;
    Q_psi_psihat(:,2) = p_w12_x_psi;
    Q_psi_psihat(:,3) = p_w2_x_psi;
    
    
    %% Maximitation step
    alpha_old_1 = alpha_1;
    alpha_old_2 = alpha_2;
    
    delta_1_old = delta_1;
    delta_2_old = delta_2;
    
    nu_old = nu;
    delta_old = delta;
    
    alpha_1 = (1/N)*sum(p_w11_x_psi);
    alpha_2 = (1/N)*sum(p_w12_x_psi);
    
    delta_1 = sqrt((sum(p_w11_x_psi.*(ro.^2)))/(2*(sum(p_w11_x_psi))));
    delta_2 = sqrt((sum(p_w12_x_psi.*(ro.^2)))/(2*(sum(p_w12_x_psi))));
    
    I_0 = besseli(0,(ro.*nu_old)./(delta_old^2));
    I_1 = besseli(1,(ro.*nu_old)./(delta_old^2));
    nu = ((sum(p_w2_x_psi.*(I_1./I_0).*ro))/(sum(p_w2_x_psi)));
    
    aux = (ro.^2) + (nu_old.^2) - (2*ro.*nu_old.*(I_1./I_0));
    delta = sqrt((sum(p_w2_x_psi.*aux))/(2*(sum(p_w2_x_psi))));
    
    
    %% Log-likelihood
    p_x_w11 = rayleigh(ro,delta_1);
    p_x_w12 = rayleigh(ro,delta_2);
    p_x_w2 = rician(ro,delta,nu);
    
    %Solve problem with low probabilities round to zero in order to avoid
    %Nan or Inf due to de fisr order modified Bessel function and
    %log-likelihood
    p_x_w11(p_x_w11 == 0) = eps;
    p_x_w12(p_x_w12 == 0) = eps;
    p_x_w2(p_x_w2 == 0) = eps;
    p_x_w2(isnan(p_x_w2)) = eps;
    p_x_w2(isinf(p_x_w2)) = eps;
    
    %joint distribution
    p_ro_Psi = (alpha_1*p_x_w11) + ...
        (alpha_2*p_x_w12) + ...
        ((1 - alpha_1 - alpha_2)*p_x_w2);
    
    log_likeliold = log_likelihood(k);
    
    k = k + 1;
    
    log_likelihood(k) = sum(log(p_ro_Psi));
    
    %     plot(log_likelihood,'r','linewidth',3)
    %     title('Log-likelihood','FontSize',13)
    %     xlabel('Iterations','FontSize',13)
    %     ylabel('$\mathcal{L}(\mathbf{ \Psi }) = \sum log(p(\mathbf{X}|\mathbf{\Psi}))$','FontSize',13,'interpreter','latex')
    %     set(gca,'FontSize',12)
    %     grid on
    %
    %     drawnow
    
end

plot(log_likelihood,'r','linewidth',3)
title('Log-likelihood','FontSize',13)
xlabel('Iterations','FontSize',13)
ylabel('$\mathcal{L}(\mathbf{ \Psi }) = \sum log(p(\mathbf{X}|\mathbf{\Psi}))$','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
grid on

%% BDR (Bayesian Decision Rule)
[val , index] = max(Q_psi_psihat,[],2);
clear val;

change_map = ro;
change_map(index == 1) = 0; %No change
change_map(index == 2) = 0; %No change
change_map(index == 3) = 1; %Change
change_map = reshape(change_map,m_c,n_c);

figure, imshow(change_map), colorbar

[counts ,centers] = hist(ro,nbins);
counts_n = counts./max(counts);
figure, bar(centers,counts_n)
hold on
plot(ro,p_ro_Psi./max(p_ro_Psi),'rx')
title('Joint distribution and Histogram of $\rho$','FontSize',13,'interpreter','latex')
xlabel('$\rho$','FontSize',13,'interpreter','latex')
ylabel('Normalized probabilities','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
grid on
