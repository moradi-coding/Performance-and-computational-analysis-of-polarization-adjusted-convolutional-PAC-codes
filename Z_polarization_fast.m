function ZW = Z_polarization_fast(N,sigma)
% SNR_dB = SNR_vec(i); % 2.5dB
% ebn0 = 10^(SNR_dB/10); % R = .5;
% sigma = 1/sqrt(2*R*ebn0); % sigma = 0.7499;

fc = @(t) phi_inv(1-(1-phi_fun(t))^2); % Check node(left branch) operation
fv = @(t) 2*t; % Bit node(right branch) operation

n = log2(N);
m = zeros(1,2*N-1); % all the tree nodes
m(1) = 2/sigma^2; % Ititialize the tree root mean
% In the tree, the left branch is check operation(fc) and right branch is
% node operation(fv)
for d = 1:n
    u = 0;
    for i = 0:2^d-1
       if u == 2^d
           break; % nodes in the depth d finished          
       end
       if mod(i+1,2) == 1 % left branch operation
          m(2^d+i) = fc(m(2^(d-1)+floor(i/2))) ;
       else % right branch operation
           m(2^d+i) = fv(m(2^(d-1)+floor(i/2)));
       end
       u = u+1;
    end
end

ZW = exp(-m(N:2*N-1)/4); % Z(W) = e^(-1/2sigma^2), m = 2/sigma^2;

end

function phi = phi_fun(t)
% The phi function in the Gaussian approximation
load('phiF_vec.mat');
if t == 0
    phi = 1;
elseif t < 50
    [~,I] = min(abs(A-t));
    phi = phi_vec(I) ;
%     syms z
%     F = tanh(z/2).*exp(-(z-t).^2./(4*t));
%     phi = 1 - 1./(sqrt(4*pi*t)).*double(int(F,z,-100,100));
else
    phi = 1.1795e-06;
end

end

function t = phi_inv(y)
load('phi_vector2.mat'); % contains x = .01:.01:45; and corresponding phi values
if y == 1
    t = 0;    
else
    [~,I] = min(abs(phi_vector-y));
    t = x(I) ;
end
end

% The way I calculated phi_vector.mat 
% x = .01:.01:45;
% syms z
% F = tanh(z/2).*exp(-(z-x).^2./(4*x));
% Integ = int(F,z,-100,100);
% Integ = double(Integ);
% phi_vector = 1 - 1./(sqrt(4*pi*x)).*Integ;
% semilogy(x,phi_vector,'-b');
% % Q_1 = exp(-0.4527*x^0.86 + 0.0218); % 0<t<10 % GA Approximations
% % Q2 = sqrt(pi./x).*exp((-x./4).*(1-20./(7*x))); t>10 % GA Approximations




% 
% % The way I calculated phi_vector for inverse.mat 
% x = .01:.1:45;
% syms z
% F = tanh(z/2).*exp(-(z-x).^2./(4*x));
% Integ = int(F,z,-100,100);
% Integ = double(Integ);
% phi_vector = 1 - 1./(sqrt(4*pi*x)).*Integ;
% semilogy(x,phi_vector,'-b');
% % Q_1 = exp(-0.4527*x^0.86 + 0.0218); % 0<t<10 % GA Approximations
% % Q2 = sqrt(pi./x).*exp((-x./4).*(1-20./(7*x))); t>10 % GA Approximations




% 
% % the way I calculated phiF_vec.mat for forward phi
% clear; clc; close all
% A = 0.1:0.1:50;
% phi_vec = zeros(1,length(A));
% parfor i = 1:length(A)
%     phi_vec(i) = phi_fun(A(i));
%     
% end
% 
% plot(A,phi_vec)