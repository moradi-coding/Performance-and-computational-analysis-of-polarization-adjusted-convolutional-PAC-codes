function metric = calc_met_gaussian(LLR,u,b)
% Compute the Fano metric for conv coded data through polarized channel for
% BPSK modulated signal with the following mapping
% 0 --> +1
% 1 --> -1
% LLR = soft output of SC decoder
% u = hypothesis
% b = bias
Lj = exp(LLR);
% if u == 0
%     metric = log2(2/(1+1/Lj)) - b;
% else
%     metric = log2(2/(1+Lj)) - b;
% end

if u == 0
    metric = 1 -  log2(1+1/Lj) - b;
else
    metric = 1 - log2(1+Lj) - b;
end

