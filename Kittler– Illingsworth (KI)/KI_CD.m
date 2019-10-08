clear all
close all
clc

%% Initialization

%________________ lake data ___________
%load('mulargia_95_96.mat')
load('fire_2013.mat')

I = lake';
Aspot = lakef';
clear lake lakef;

%________________ contest data ___________
%load('data_contest.mat')


%_________________Magnitude of difference _______________


ro = sqrt((I(:) - Aspot(:)).^2);

T = kittler(ro);

idx_w1 = ro <= T;
idx_w2 = ro > T;
W1 = ro(idx_w1);
W2 = ro(idx_w2);

change_map = ro;
change_map(idx_w1) = 0;
change_map(idx_w2) = 1;
[m_c n_c] = size(Aspot);
figure, imshow(reshape(ro,m_c,n_c)), colorbar
figure, imshow(reshape(change_map,m_c,n_c)), colorbar