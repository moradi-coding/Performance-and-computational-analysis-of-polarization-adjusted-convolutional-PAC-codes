function x = polarencode(u)
N = length(u);
n = log2(N);
x = u;
m =1; % #of bits combined
for d = n-1:-1:0
    for i = 1:2*m:N
        a = x(i:i+m-1); % first part
        b = x(i+m:i+2*m-1); % second part
        x(i:i+2*m-1) = [mod(a+b,2) b]; % Combining
    end
    m = m*2;
end
end