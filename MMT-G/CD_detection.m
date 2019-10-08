clear all
close all
clc

%% Data read and parameters

% sampling
addpath(genpath([pwd , '\RR_ImageSampling\']));
% MI
addpath(genpath([pwd , '\minf\']));
%Data
addpath(genpath([pwd , '\Data\']));

%addpath(genpath('J:\Google Drive\altmany-export_fig-4703a84\'));


load('mulargia_95_96.mat')
%load('fire_2013.mat')


% ------------ for lake & fire data----------
Aspot = (lakef');
maxA = max(max(Aspot(:,:,1)));
Aspot = (Aspot./maxA);

I = (lake');
maxA = max(max(I(:,:,1)));
I = (I./maxA);

%%

[m_c n_c] = size(Aspot);

modalities = 2;

%Numeber of samples
n = 92;

Xl_AA = cell(modalities,1);
Xl = cell(modalities,1);

%% Extract random pixels from Images
%-----------if you prefer to get the locations manually--------------------

%imshow(Aspot(:,:,1),[]);

%[xi,yi] = getpts;
%locations = round([xi yi]);

% Manual
%locations = sub2ind([m_c  n_c],locations(:,2),locations(:,1)); %manual

% --------------Get n random locations (Automatically)--------------------
[ locations, v_j ] = Sampling_Grid(I(:,:,1), n, true );
%[ locations, v_j ] = Sampling_Jittered(I(:,:,1), n, true);
%[locations, v_j] = Sampling_Uniform(I(:,:,1), n, true );
%[ locations, v_bc ] = Sampling_BestCandidate(I(:,:,1), n, ceil(n/4), true );

%Autom
locations = sub2ind([m_c  n_c],locations(2,:),locations(1,:));


n = length(locations);

% Get values at those locations:

Xl{1} = (I(:,:,1)); %Ir
Xl{2} = (Aspot(:,:,1)); %Ig

Xl_AA{1} = Xl{1}(locations); %Ir
Xl_AA{2} = Xl{2}(locations); %Ig

% local graph
wl = cell(modalities,1);

complement = setdiff(1:(m_c*n_c), locations);


labels = {'NIR Bf','NIR Af'};


%------another metric taking into account the index of the pixel

% indexdisa = pdist2(locations',locations');
% max_distance = max(nonzeros(indexdisa));
% indexdisa = indexdisa./max_distance;
% indexistda = std(indexdisa(:));
% % %
% indexdisb = pdist2((complement)',locations');
% max_distance = max(nonzeros(indexdisb));
% indexdisb = indexdisb./max_distance;
% indexistdb = std(indexdisb(:));
% 
% metrica = @(x,kernelstda) exp(-x ./ (2*(kernelstda ^ 2))).*...
%  exp(-indexdisa ./ (2*(indexistda ^ 2)));
% %
% metricb = @(x,kernelstda) exp(-x ./ (2*(kernelstda ^ 2))).*...
% exp(-indexdisb ./ (2*(indexistdb ^ 2)));

%-------------Kernel heat-------------------------
metric = @(x,kernelstda) exp(-(x.^2) ./ ((kernelstda ^ 2)));


for i = 1 : modalities
    
    %% L2 norm
    %[ind, distl{i}] = knnsearch(Xl{1}(:), Xl{1}(:),...
    %     %'K', 10, 'Distance', 'euclidean');
    %distlAA = squareform(pdist(Xl_AA{i}(:)));
    distlAA = pdist2(Xl_AA{i}',Xl_AA{i}','euclidean');
    distlAB = pdist2(Xl{i}(complement)',Xl_AA{i}','euclidean').^3;
    
    
    %% L1 norm
%                 x=repmat(Xl_AA{i}',[size(Xl_AA{i}',1),1]);
%                 y=repmat(Xl_AA{i},[size(Xl_AA{i},2),1]);
%                 y = y(:);
%                 distlAA = reshape(sum(abs(x-y),2).^2,n,n);
%                 x=repmat(Xl_AA{i},[size(Xl{i}(complement),2),1]);
%                 x = x(:);
%                 y=repmat(Xl{i}(complement)',[size(Xl_AA{i}',1),1]);
%                 distlAB = reshape(sum(abs(x-y),2).^2,length(complement),n);
    %
    %% Normalization and Degree
    
    D1 = distlAA*ones(n,1) + distlAB'*ones(length(complement),1);%diag(sum(distlAA, 2));
    distlAA = distlAA./repmat((D1),1,n);
    D2 = distlAB*ones(n,1) + (distlAB*pinv(distlAA))*(distlAB'*ones(length(complement),1));%sum(distlAB, 2);
    distlAB = distlAB./repmat((D2),1,n);
    
    
    clear x y D1 D2
    
    %% Metric computation
    
    %auxstd = [distlAA ; distlAB];
    %mean data
    kernelstd = mean(distlAB(:));
    %silverman
    %kernelstd = 2.07*std(distlAB(:))/(n^(1/5));
    %clear auxstd
    
    distlAA = metric(distlAA,kernelstd);
    distlAB = metric(distlAB,kernelstd);
    
    
    wl{i} = [distlAA ;distlAB];
    
    clear distlAA distlAB distlAAsmooth
    
    figure
    imagesc(wl{i}), colorbar;
    title(['Weights of ' labels{i} ' band'])
    set(gca,'FontSize',12)
    
    drawnow;
end
clear Xl Xl_AA
%% Multimodal Weights

n = length(locations);
WL = min(cat(3,wl{1} , wl{2}),[],3);
figure
imagesc(WL), colorbar
title(['Multimodal Weights'])
set(gca,'FontSize',12)

% One shot method of Nystrom
W_AA = WL(1:n,1:end);
W_BA = WL(n+1:end,1:end);


%% Nystrom aproximation

W_AA_sqrtinv = pinv(sqrtm(W_AA));%W_AA^-0.5;
%
S = W_AA + (W_AA_sqrtinv*(W_BA'*W_BA)*W_AA_sqrtinv);

[U_s,D_s] = eig(S);

Uhat_W = [ W_AA ; (W_AA_sqrtinv*W_BA')']*(U_s*pinv(sqrtm(D_s)));%%%[U_s; W_BA*U_s*(pinv(D_s))];%[W_AA;W_BA*W_AA_sqrtinv]*(U_s*pinv(sqrtm(D_s)));%%%

clear S W_AA_sqrtinv W_AA W_BA

sel = n - 2;
figure, imshow(reshape(Uhat_W(:,sel)*sqrt(D_s(sel,sel)),m_c,n_c),[]),colorbar%, colormap hot
set(gca,'FontSize',12)

h = figure; imshow(I(:,:,1),[])%, colorbar
title('Before')
set(gca,'FontSize',12)
%export_fig(strcat(pwd,'\Figures\bf_fire'),'-pdf','-transparent',h)


h = figure; imshow((Aspot(:,:,1)),[])%, colorbar
title('After')
set(gca,'FontSize',12)
%export_fig(strcat(pwd,'\Figures\af_fire'),'-pdf','-transparent',h)


%E = zeros(1,n);
MI = zeros(1,n);
prior = ((I(:,:,1) - Aspot(:,:,1))./(I(:,:,1) + Aspot(:,:,1)));
prior = imbinarize(prior);
figure, imshow((prior(:,:,1)),[]), colorbar
title('Prior')
set(gca,'FontSize',12)


for i = 1 : n
    Iaux = Uhat_W(:,i)*sqrt(D_s(i,i));
    
    A = Iaux(1:n);
    AB = Iaux(n+1:end);
    Iaux(locations) = A;
    Iaux(complement) = AB;
    Iaux = ((reshape(Iaux,m_c,n_c)));
    
    Iaux = imbinarize(Iaux);
    
%     figure, imshow(Iaux,[]),colorbar%, colormap hot
%     title(['Eigenvector sample ' , num2str(i) ])
%     set(gca,'FontSize',12)
    
    MI(i) =  mi(prior,Iaux);
    clear Iaux;
end

figure, plot(MI,'linewidth',2)
title('MI of eigenvectors')
xlabel('$u_{i} \sqrt{d_{i}}$','FontSize',13,'interpreter','latex')
ylabel('$MI(I_{u_{i}},I_{Prior})$','FontSize',13,'interpreter','latex')
set(gca,'FontSize',12)

i = find(MI == max(MI),1,'first');

hold on

plot(i,MI(i),'ro','MarkerSize',12)

Iaux = Uhat_W(:,i)*sqrt(D_s(i,i));

A = Iaux(1:n);
AB = Iaux(n+1:end);
Iaux(locations) = A;
Iaux(complement) = AB;
Iaux = ((reshape(Iaux,m_c,n_c)));

Iaux = imbinarize(Iaux);

figure, imshow(Iaux,[])
title(['Eigenvector sample ' , num2str(i) ])
set(gca,'FontSize',12)
