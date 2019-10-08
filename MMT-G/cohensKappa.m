function [MA, FA, P, R, kappa, OE] = cohensKappa(y, yhat)
    C = confusionmat(y, yhat); % compute confusion matrix
    
    %% computation of MA FA precision, recall and OE
    total_cd_pixels = length(find(y));
    total_ncd_pixels = length(y) - total_cd_pixels;
    
    MA = (C(2,1)/total_cd_pixels)*100;
    FA = (C(1,2)/total_ncd_pixels)*100;
    
    P = C(2,2)/(C(2,2) + C(1,2));%1 - FA/(total_cd_pixels - MA + FA);
    R = C(2,2)/(C(2,2) + C(2,1));%1 - (MA/total_cd_pixels);
    
    OE = ((C(1,2) + C(2,1))/sum(C(:)))*100;
    
    %% computation of kappa
    n = sum(C(:)); % get total N
    C = C./n; % Convert confusion matrix counts to proportion of n
    r = sum(C,2); % row sum
    s = sum(C); % column sum
    expected = r*s; % expected proportion for random agree
    po = sum(diag(C)); % Observed proportion correct
    pe = sum(diag(expected)); % Proportion correct expected
    kappa = (po-pe)/(1-pe); % Cohen's kappa
end