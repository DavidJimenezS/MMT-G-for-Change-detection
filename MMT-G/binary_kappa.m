clear all
close all
clc

addpath(genpath([pwd , '\Data\']));
addpath(genpath([pwd , '\Results\']));
addpath(genpath('J:\Google Drive\altmany-export_fig-4703a84\'));


warning off;
%% Correcting locations of the samples pixels 
%-----------------Data-------------------------
load('lake_95_96_cub.mat'), load('gt_lake.mat'), load('rR_lake.mat'), load('rrR_lake.mat')
load('KI_lake.mat')

Iaux3_nofilt = (U(:,1))*sqrt(Eval);

n = length(locations);
A = Iaux3_nofilt(1:n);
AB = Iaux3_nofilt(n+1:end);
Iaux3_nofilt(locations) = A;
Iaux3_nofilt(complement) = AB;
Iaux3_nofilt = ((reshape(Iaux3_nofilt,m_c,n_c)));

clear U Eval m_c n_c locations complement A AB


Ibw1_lake = imbinarize(Iaux3_nofilt);
figure, imshow(Ibw1_lake,[])
title('Change map using MMG','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)


rR_bw_lake = imbinarize(change_map);
figure, imshow(rR_bw_lake,[])
title('Change map using rR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

rrR_bw_lake = imbinarize(rrR_change_map);

figure, imshow(rrR_bw_lake,[])
title('Change map using rrR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

KI_bw_lake = imbinarize(KI_lake_change_map);

figure, imshow(KI_bw_lake,[])
title('Change map using KI','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)


gt_bw = imbinarize(gt_lake);
h = figure; imshow(gt_bw,[])
title('Reference Change map','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
%export_fig(strcat(pwd,'\Figures\ref_lake'),'-pdf','-transparent',h)

MA = zeros(4,1); FA = MA; Precision = MA; recall = MA; kappa = MA; OE = MA;
methods = {'rR-EM';'rrR-EM';'KI';'MMG'};

disp('Results for lake Mulargia flood event')

[MA(1), FA(1), Precision(1), recall(1), kappa(1), OE(1)] = cohensKappa(gt_bw(:),rR_bw_lake(:));

[MA(2), FA(2), Precision(2), recall(2), kappa(2), OE(2)] = cohensKappa(gt_bw(:),rrR_bw_lake(:));

[MA(3), FA(3), Precision(3), recall(3), kappa(3), OE(3)] = cohensKappa(gt_bw(:),KI_bw_lake(:));

[MA(4), FA(4), Precision(4), recall(4), kappa(4), OE(4)] = cohensKappa(gt_bw(:),Ibw1_lake(:));

T1 = table(MA,FA,Precision,recall,kappa,OE,'RowNames',methods)


%% ---------------------------------- ----------------

%-----------------Data set 2---------------------


%% Correcting locations of the samples pixels
%-----------------Data-------------------------
load('fire_cub.mat'), load('ref_fire.mat'),  load('rR_fire.mat'), load('rrR_fire.mat')

load('KI_fire.mat')

Iaux3_nofilt = (U(:,1))*sqrt(Eval);

n = length(locations);
A = Iaux3_nofilt(1:n);
AB = Iaux3_nofilt(n+1:end);
Iaux3_nofilt(locations) = A;
Iaux3_nofilt(complement) = AB;
Iaux3_nofilt = ((reshape(Iaux3_nofilt,m_c,n_c)));

clear U Eval m_c n_c locations complement A AB


%% Fire results

Ibw1_fire = imbinarize(Iaux3_nofilt);
figure, imshow(Ibw1_fire,[])
title('Change map using MMG 3','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

rR_bw_fire = imbinarize(rR_fire_map);
figure, imshow(rR_bw_fire,[])
title('Change map using rR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

rrR_bw_fire = imbinarize(rrR_fire_map);
figure, imshow(rrR_bw_fire,[])
title('Change map using rrR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

KI_bw_fire = imbinarize(KI_fire_change_map);
figure, imshow(KI_bw_fire,[])
title('Change map using KI','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

gtf_bw = imbinarize(gt_fire);
h = figure; imshow(gtf_bw,[])
title('Reference Change map','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
%export_fig(strcat(pwd,'\Figures\ref_fire'),'-pdf','-transparent',h)


MA = zeros(4,1); FA = MA; Precision = MA; recall = MA; kappa = MA; OE = MA;
%methods = {'MMG';'rR-EM';'inverse rR-EM'};
methods = {'rR-EM';'rrR-EM';'KI';'MMG'};


%I3_bw = imbinarize(Iaux3);

[MA(1), FA(1), Precision(1), recall(1), kappa(1), OE(1)] = cohensKappa(gtf_bw(:),rR_bw_fire(:));

[MA(2), FA(2), Precision(2), recall(2), kappa(2), OE(2)] = cohensKappa(gtf_bw(:),rrR_bw_fire(:));

[MA(3), FA(3), Precision(3), recall(3), kappa(3), OE(3)] = cohensKappa(gtf_bw(:),KI_bw_fire(:));

[MA(4), FA(4), Precision(4), recall(4), kappa(4), OE(4)] = cohensKappa(gtf_bw(:),Ibw1_fire(:));

disp('Results for fire event')

T2 = table(MA,FA,Precision,recall,kappa,OE,'RowNames',methods)


%% Show map with respect to MA, FA and correct detection

%---------lake------------------------------------------------------
%rR
E_lake_rR = gt_bw - rR_bw_lake;
E_lake_rR(E_lake_rR == -1) = 3; %FA
E_lake_rR(E_lake_rR == 1) = 2; %MA
aux_1 = gt_bw + rR_bw_lake;
E_lake_rR(aux_1 == 2)= 1;

rR_label_lake = label2rgb(E_lake_rR,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(rR_label_lake)
title('Error map using rR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\rR-EM_Emap_lake'),'-pdf','-transparent',h)

%rrR
E_lake_rrR = gt_bw - rrR_bw_lake;
E_lake_rrR(E_lake_rrR == -1) = 3; %FA
E_lake_rrR(E_lake_rrR == 1) = 2; %MA
aux_1 = gt_bw + rrR_bw_lake;
E_lake_rrR(aux_1 == 2)= 1;

rrR_label_lake = label2rgb(E_lake_rrR,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(rrR_label_lake)
title('Error map using rrR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\rrR-EM_Emap_lake'),'-pdf','-transparent',h)


%KI
E_lake_KI = gt_bw - KI_bw_lake;
E_lake_KI(E_lake_KI == -1) = 3; %FA
E_lake_KI(E_lake_KI == 1) = 2; %MA
aux_1 = gt_bw + KI_bw_lake;
E_lake_KI(aux_1 == 2)= 1;

KI_label_lake = label2rgb(E_lake_KI,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(KI_label_lake)
title('Error map using KI','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\KI_Emap_lake'),'-pdf','-transparent',h)


%MMG
E_lake_MMG = gt_bw - Ibw1_lake;
E_lake_MMG(E_lake_MMG == -1) = 3; %FA
E_lake_MMG(E_lake_MMG == 1) = 2; %MA
aux_1 = gt_bw + Ibw1_lake;
E_lake_MMG(aux_1 == 2)= 1;

MMG_label_lake = label2rgb(E_lake_MMG,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(MMG_label_lake)
title('Error map using MMG','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\MMG_Emap_lake'),'-pdf','-transparent',h)

%-------------------Fire-----------------
%rR
E_fire_rR = gtf_bw - rR_bw_fire;
E_fire_rR(E_fire_rR == -1) = 3; %FA
E_fire_rR(E_fire_rR == 1) = 2; %MA
aux_1 = gtf_bw + rR_bw_fire;
E_fire_rR(aux_1 == 2)= 1;

rR_label_fire = label2rgb(E_fire_rR,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(rR_label_fire)
title('Error map using rR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\rR-EM_Emap'),'-pdf','-transparent',h)


%rrR
E_fire_rrR = gtf_bw - rrR_bw_fire;
E_fire_rrR(E_fire_rrR == -1) = 3; %FA
E_fire_rrR(E_fire_rrR == 1) = 2; %MA
aux_1 = gtf_bw + rrR_bw_fire;
E_fire_rrR(aux_1 == 2)= 1;

rrR_label_fire = label2rgb(E_fire_rrR,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(rrR_label_fire)
title('Error map using rrR-EM','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\rrR-EM_Emap'),'-pdf','-transparent',h)


%KI
E_fire_KI = gtf_bw - KI_bw_fire;
E_fire_KI(E_fire_KI == -1) = 3; %FA
E_fire_KI(E_fire_KI == 1) = 2; %MA
aux_1 = gtf_bw + KI_bw_fire;
E_fire_KI(aux_1 == 2)= 1;

KI_label_fire = label2rgb(E_fire_KI,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(KI_label_fire)
title('Error map using KI','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\KI_Emap'),'-pdf','-transparent',h)


%MMG
E_fire_MMG = gtf_bw - Ibw1_fire;
E_fire_MMG(E_fire_MMG == -1) = 3; %FA
E_fire_MMG(E_fire_MMG == 1) = 2; %MA
aux_1 = gtf_bw + Ibw1_fire;
E_fire_MMG(aux_1 == 2)= 1;

MMG_label_fire = label2rgb(E_fire_MMG,[0 1 0;0 0 1; 1 0 0],[1 1 1]);
h = figure; imshow(MMG_label_fire)
title('Error map using MMG','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)
export_fig(strcat(pwd,'\Figures\MMG_Emap'),'-pdf','-transparent',h)
