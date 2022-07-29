function RP = RM_Polar_Profile(N, K, SNR, type)
% type 1 -------> GA consttruction, see [2]
% type 2 -------> RMpolar construction see [1].
% type 3 -------> Tse RMpolar construction
% type 4 -------> RM profile, see [4, 2]
% N = 1024; K = 512; SNR = 2; 
R = K/N;

EbNo_X = 10^(SNR/10); 
sigma_X = 1/sqrt(2*R*EbNo_X);
ZW_X = Z_polarization_fast(N,sigma_X);

if type == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRP rate profile (K most reliable channels)
    RP = false(1,N);    
    [~,I] = mink(ZW_X, K);
    RP(I) = true;
%     sum(RP)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Klow1 = 0;
    n = log2(N);
    i = n;
    while Klow1 < K    
        if Klow1 + nchoosek(n,i) < K
            Klow1 = Klow1 + nchoosek(n,i);
        else
            break;
        end
        i = i-1;
    end
    
    [RP, w] = RM_profile(Klow1,N);
    idx_w = find(w==i); % rows with weight i
           
    E0_X = log2(2./(1+ZW_X));
    
    X = zeros(1,N);
    X(idx_w) = E0_X(idx_w);
    [~,I] = sort(X,'descend');  
    I(K-Klow1+1:end) = []; % most reliable rows with weight i
    RP(I) = true;
%     sum(RP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Klow2 > K and we choose the most reliable ones from it
    Klow2 = 0;
    n = log2(N);
    i = n;
    while Klow2 <= K    
        Klow2 = Klow2 + nchoosek(n,i);
        i = i-1;
    end
    [RP, ~] = RM_profile(Klow2,N);
        
    E0_X = log2(2./(1+ZW_X));    
    
    E = 1e6*E0_X;
    RP_X = RP.*E;
    [~,I] = maxk(RP_X, K);
    
    RP = false(1,N);
    RP(I) = true;
%     sum(RP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif type == 4
    [RP, ~] = RM_profile(K,N);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, w] = RM_profile(K,N)
w = sum(dec2bin(0:N-1)-'0',2);
[~, index] = sort(w,'descend');
kindex = index(1:K);
P = false(1,N);
P(kindex) = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









