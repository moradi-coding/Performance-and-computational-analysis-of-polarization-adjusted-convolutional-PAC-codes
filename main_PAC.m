clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for "parfor" and CoreNum = 1 is like using simple "for"
CoreNum = 4;
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(CoreNum);
else
    disp('matlab pool already started');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = rng(100);       % Seed the RNG for repeatability

% delta is a Fano decoder parameter whose effect on decoding PAC codes is 
% described in [2]. If the complexity is not finite, you may use 2. 
delta = 2; 
% You may use a finite #cycles in practice or large N 
maxcycles = 1e15;

EbNo_vec = 1:0.5:3; % we have loop for different SNR values


% SNR_Cons is the construction SNR for Gaussian app, or RM-polar construction. 
% We use RM construction in [2,3 ...]. 
% An MC-based constrction is introduced in [5] based on the channel
% polarization. You may use or introduce any Cons that works well for PAC codes.
K = 64; N = 128; R = K/N; Rc = R; SNR_Cons = 1;

% rate profile
type = 4; % 1--> GA, 3--> Tse_RMpolar, 4 --> RM
RP = RM_Polar_Profile(N, K, SNR_Cons, type);
sum(RP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convolutional code that we use matrix T in our references.
% We designed c(x) = 3211 as explained in [4, Ch6.2 1]. You may use your poly.
% You may implement conv encoding in a better way for longer poly
poly = 3211;
polyb = dec2bin(base2dec(num2str(poly), 8)) - '0'; % 1 011 011
constraint = length(polyb); % 7
C = constraint - 1;
tr = poly2trellis(constraint,poly); % #of states is 2^C
outputs = tr.outputs; 
nextStates = tr.nextStates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FER = zeros(1,length(EbNo_vec)); % Frame error rate for each SNR
FE = zeros(1,length(EbNo_vec));  % #of frames in errors
V = zeros(1,length(EbNo_vec));  % #of visits for every SNR
Nruns = zeros(1,length(EbNo_vec)); % #of actual runs

maxRun = 1e9;

ClustNum = [500 500 1e3 5e3*(ones(1,length(EbNo_vec)-3))]; % adjust this based on N and K
maxFE = 100; % maximum number of frame errors, you may use 300 for smoother curves

fprintf('-------------------------------------\n');
for EbNo_count = 1:length(EbNo_vec)
    tic;
    EbNo_dB = EbNo_vec(EbNo_count);
    EbNo = 10^(EbNo_dB/10);
    sigma = 1/sqrt(2*Rc*EbNo);
    
% We use bit-channel cutoff rate E0 as the bias value for the Fano decoder.
% Refere to [2] for its implementation and its affect on the complexity. In
% [4] we use fixed bias which is different for different R. In hardware you
% may 1 bit quantize the vector E0 as E0(E0>=.5)=1 & E0(E0<.5)=0.
    ZW = Z_polarization_fast(N,sigma);
    E0 = log2(2./(1+ZW));
%     IW = I_polarization(N,sigma); insted of E0 gives better FER for low
%     SNR values (higher complexity), see [2].
    
    Nblkerrs = 0;
    Visit = 0;
%     updateDisp(SNR_vec,percent,FE,FER,Nruns,i);
    fprintf('[%02d:%02d] Starting! SNR = %.2f\n',0,0, EbNo_dB);
    for i = 1:maxRun/ClustNum(EbNo_count)
        
        parfor j = 1:ClustNum(EbNo_count)
            %Generate random message
            msg = randi([0 1],1,K);
            code = zeros(1,N);
            code(RP) = msg; % rate profiling step
            % Convolutional encoding
            u = convenc(code,tr); % convolutional encoder step.
            % Polar transformation
            x = polarencode(u);
            % BPSK modulation
            modulated = 2*x - 1;
            % AWGN
            r = modulated + randn(1,N)*sigma;
            % Decoding
            [~,~,decoded,VPC] = ...
                pac_decode2(r,delta,maxcycles,poly,RP,sigma,E0);
            % You may count visits only for forward move, or forward\backward node
            % visits or f/g-like visits [2].
            Visit = Visit + VPC;
            
            Nblkerrs = Nblkerrs + any(decoded(RP)~=msg');
        end
        if Nblkerrs >= maxFE
            break;
        end
        t = toc;
        elapsed_m = t/60;
        elapsed_h = floor(elapsed_m/60);
        elapsed_m = floor(elapsed_m - elapsed_h*60);
        fprintf(2,'[%02d:%02d] EbNo = %.2f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs);
    end
    t = toc;
    elapsed_m = t/60;
    elapsed_h = floor(elapsed_m/60);
    elapsed_m = floor(elapsed_m - elapsed_h*60);
    fprintf(2,'[%02d:%02d] SNR = %.1f, Frame = %d, FE = %d\n',elapsed_h,elapsed_m,EbNo_dB,i*ClustNum(EbNo_count),Nblkerrs);
    
    temp = (i*ClustNum(EbNo_count));
    Nruns(EbNo_count) = temp;
    FER(EbNo_count) = Nblkerrs/temp;
    FE(EbNo_count) = Nblkerrs;    
    V(EbNo_count) = Visit;
    fprintf('-------------------------------------\n');
end

rng(s);     % Restore RNG

% Creating sim_state
res.trials = Nruns;
res.frame_errors = FE;
res.FER = FE./Nruns;
res.avg_visit = V./Nruns;

res.SNR = EbNo_vec;
res.rate_profile = type;
res.K = K;
res.delta = delta;

res.max_cycles_allowed = maxcycles;
res.max_trials = maxRun;
res.max_frame_errors = maxFE;


filename = sprintf('outputs/pac_K%d_N%d_SNR%0.1f_delta%0.1f_poly%d.mat',K,N,EbNo_vec(1),delta,poly);
save(filename,'res');

