% vectorized Go back
function [metric,cycles,data,VPC] = ...
    pac_decode2(r,delta,maxcycles,poly,RP,sigma,E0)
% For detailed explenation of PAC Fano decoder refere to [Ch2. 1] and [2].
% For low memery stack/heap decoder see [Ch3.5 1].
% For fast SCL (avoiding sorting) (idea from Fano decoder) see [3 (or Ch5 1)].
% For implementing stack with small stack size see [7].
% For improved (fast) stack decoder based on polarized channels [3 (or Ch5.2 1)];

% The Fano algorithm for decoding of Rate one convolutional code
% wmc: With Metric Calculator
% r = [r(0), r(1), ... r(L-1)]
% bias = fano metric bias value
% delta = Threshold increment size
% snr = current snr in decimal (not dB)
% maxcycles = Decoding timeout in cycles per bit
% poly: generator polynomials in octal
% RP = rate profile: a 1xL (1xK) vector forcing the decoder to chose a particular path
% ------ 0:chose path 0(frozen bit), 1: chose best path

% visits = zeros(1,length(RP));
VPC = 0;
polyb = dec2bin(base2dec(num2str(poly), 8))-'0';
polyL = polyb == 1;         % Logical array
c = length(polyb) - 1;      % Constraint Length


N = size(r,2);                          % number of branches
i = 1;                                  % node index
T = 0;                                  % threshold
M = zeros(N+2,1);                       % list of metrics
M(1) = -inf;                            % metric backward from root
M(2) = 0;                               % metric at root

t = zeros(N,1);                         % branch being tested at each node
inpseq = zeros(N,1);                    % input sequence to current node
parseq = zeros(N,1);                    % convolutional parity sequence to current node
tm = zeros(2,1);                        % sorted metrics for current hypotheses
tb = zeros(2,1);                        % current hypotheses sorted by corresponding metric
parity = zeros(2,1);                    % current parity sorted by corresponding metric

reg = zeros(1,c + N);                   % shift register (Conv encoder replica in the decoder

maxcycles = maxcycles * N;              % bitwise maximum cycle
dobreakA = 0;

% Metric calculator parameters
n = log2(N);
llr = zeros(1,N);       % output beliefs
L = zeros(n+1,N);       % node beliefs
ucap = zeros(n+1,N);    % decisions
ns = zeros(1,2*N-1);    % node state vector

Tanf = @(x) abs(log(tanh(abs(x)/2)));
f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*Tanf(Tanf(abs(a))+Tanf(abs(b)));
% f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b));
g = @(a,b,c) b+(1-2*c).*a;
L(1,:) = -2*(1/sigma^2)*r;      % belief of root
node = 0; depth = 0;    % start at root
% done = 0;             % decoder has finished

cycles = 0;

while cycles <= maxcycles % loop A
    t(i) = 1; % branch with the best metric is chosen (t(i) = 2, second best)
    if(dobreakA)
        break;
    end
    while(1) % loop B
        cycles = cycles + 1;
        if i - 1 < node  % i is starting from one, -1 make it start from 0
            % we node to go back, we should update the states like we have not
            % gone more than node i-1, we sould erase future from node i-1
            
            ucap(n+1,:) = parseq;
            %%%%%%%%%%%%%%%%%%%%%%Erasing Future%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % when we go up, states should be zero as the case we have not
            % touched these nodes
            
            for I = (node+1):-1:i+1 % going back in the leaf node
                npos = (2^depth-1) + I;
                
% % %                 %%%%%%%%%%%%%%%%%%%%%%%%%vectorize going back%%%%%%%%%%%%%
% % %                 v = 2.^(-(0:n));
% % %                 Nozero = 1- logical((v*npos*2)-floor(v*npos*2) );
% % %                 Future = nonzeros ( floor( (v*npos).*Nozero ) )';
% % %                 ns(Future(1:end-1)) =  0;
% % %                 ns(Future(end))= 1;
% % %                 %%%%%%%%%%%%%%%%%%%%%%%%%vectorize going back rnd%%%%%%%%%
                
                ns(npos) = 0; % always we start with going up
                tempU = npos;
                while (mod(tempU,2) == 0) % condition which we should go up more
                    tempU = tempU/2;
                    ns(tempU) = 0;
%                     ns(tempU-1) = 2; % the nodes we went right ( we have tendency to
                    %go right [every node we had gone up, in its left node we go right]),
                    % in goback we had this step, but actually we dont need it.    
                end
                tempL = floor(tempU/2); %we can go up to a limit, then we start to go right.
                ns(tempL) = 1;
            end
            depth = n; % we are in the leaf
            node = (i-1); % at the end node position is updated to the place we need the llr
            
            llr(node+1) = L(n+1,node+1); % llr starts from 1
            %%%%%%%%%%%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            stop = i; % want to depolarize upto leaf node i-1 (nodes start from zero)
            ucap(n+1,:) = parseq;
            %%%%%%%%%%%%%%%%%%%%%%Depolarize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this script is used in fano_polar_wmc
            % stop must be defined
            done = 0;               % decoder has finished
            while(done == 0)
                if depth == n  % leaf or not
                    llr(node+1) = L(n+1,node+1);
                    if node == (stop-1)  % stop = i;
                        done = 1;
                    else
                        node = floor(node/2); depth = depth - 1;
                    end
                else
                    % nonleaf
                    npos = (2^depth-1) + node + 1;  % position of node in node state vector
                    if ns(npos) == 0 % step L and go to the left child
                        temp = 2^(n-depth);
                        Ln = L(depth+1,temp*node+1:temp*(node+1));  % incomin beliefs
                        a = Ln(1:temp/2); b = Ln(temp/2+1:end);
                        node = node*2; depth = depth+1; % next node: left child
                        temp = temp/2;  % incoming belief length for left child
                        L(depth+1,temp*node+1:temp*(node+1)) = f(a,b);
                        ns(npos) = 1;
                    else
                        if ns(npos) == 1 % step R and go to the right chilc
                            temp = 2^(n-depth);
                            Ln = L(depth+1,temp*node+1:temp*(node+1));  % incomin beliefs
                            a = Ln(1:temp/2); b = Ln(temp/2+1:end);
                            lnode = 2*node; ldepth = depth+1; % left child
                            ltemp = temp/2;
                            ucapn = ucap(ldepth+1,ltemp*lnode+1:ltemp*(lnode+1));
                            node = node*2+1; depth = depth+1; % next node: right child
                            temp = temp/2;  % incoming belief length for right child
                            L(depth+1,temp*node+1:temp*(node+1)) = g(a,b,ucapn);
                            ns(npos) = 2;
                        else % Step U and go to the parent
                            temp = 2^(n-depth);
                            lnode = 2*node; rnode = 2*node+1; cdepth = depth+1; % left and right child
                            ctemp = temp/2;
                            ucapl = ucap(cdepth+1,ctemp*lnode+1:ctemp*(lnode+1));
                            ucapr = ucap(cdepth+1,ctemp*rnode+1:ctemp*(rnode+1));
                            ucap(depth+1,temp*node+1:temp*(node+1)) = [mod(ucapl+ucapr,2) ucapr];
                            node = floor(node/2); depth = depth - 1;
                            ns(npos) = 3;
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%End%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        ri = llr(i);
        state = reg(1:c);% state at current node, by moving forward/backward state will change
        
        if ~RP(i) % Frozen bit
%             visits(i) = visits(i) + 1;
            u0 = genparity_Rone_logical(0,state,polyL); % conv encovder
            m = calc_met_gaussian(ri,u0,E0(i));
            M(i+2) = M(i+1) + m;
            inpseq(i) = 0;
            parseq(i) = u0; % we show conv output with u
        else      % information bit
%             visits(i) = visits(i) + 1;
            u0 = genparity_Rone_logical(0,state,polyL);
            m0 = calc_met_gaussian(ri,u0,E0(i)); % bit metric
            
            u1 = genparity_Rone_logical(1,state,polyL);
            m1 = calc_met_gaussian(ri,u1,E0(i));
            
            if m1 > m0
                tm(1) = m1;
                tb(1) = 1;
                parity(1) = u1; % we show conv output with u
                tm(2) = m0;
                tb(2) = 0;
                parity(2) = u0;
            else
                tm(1) = m0;
                tb(1) = 0;
                parity(1) = u0;
                tm(2) = m1;
                tb(2) = 1;
                parity(2) = u1;
            end
            
            % M(i+2) = forwards node metric; Mf
            % M(i+1) = current node metric; M
            M(i+2) = M(i+1) + tm(t(i));
            inpseq(i) = tb(t(i));
            parseq(i) = parity(t(i));
        end
        
%         if M(i+1) < T + delta %if the previous node violates T+delta, currenct node has visited for the first time
%             if RP(i) % Frozen bit
%                 visits_Inf = visits_Inf + 1;   
%             end
%             visits = visits + 1;   
%         end
        
        if M(i+2) >= T % Node is acceptable (Mf >= T)
            if M(i+1) < T + delta % First time visit; tighten threshold
%                 visits(i) = visits(i) + 1;
%                     visits = visits + 1;

                while T + delta <= M(i+2)
                    T = T + delta;
                end
            end
            % Move forward
            %register shifted to right and the first reg is replaces by input
            reg = [inpseq(i) reg(1,1:end-1)];
            i = i + 1;
            VPC = VPC + 1;
            if i == N + 1 % we have reached the end of sequence
                dobreakA = 1;
            end
            break         % go to loop A
        else  % Metric does not exceed threshold
            dobreak = 0;
            while(1) % loop C
                if M(i) <  T
                    T = T - delta;
                    dobreak = 1;
                    break
                else
                    % Move back
                    i = i - 1;
                    reg = [reg(2:end) 0];%for the next input reg[end] will be replaced by reg[n-1], so no difference what we write in reg[end]
                    t(i) = t(i) + 1;
                    if ~RP(i)    % non-branching node -> look back
                        continue
                    end
                    if t(i) == 3 % worst node -> look back
                        
                        continue;       % go to loop C
                    end
                    break               % go to loop B (since dobreak = 0)
                end
            end % end of loop C
            if(dobreak)
                break;                  % break out of loop B and go to loop A
            end
        end
        
    end % end of loop B
end
data = inpseq;
% metric = M(i+1);
metric = M;

