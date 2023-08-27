rho_AB = [1/3,0,0,0,1/3,0,0,0,1/3; 
    0,0,0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0,0,0; 
    1/3,0,0,0,1/3,0,0,0,1/3; 
    0,0,0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0,0,0; 
    0,0,0,0,0,0,0,0,0; 
    1/3,0,0,0,1/3,0,0,0,1/3];   

N = 6;   % number of childrens
k = 3;    % number of outcomes in a children
dA=3;   % dimension of system A
dB=3;

c_GM = 1;

X1 = [0, c_GM, 0; c_GM, 0, 0; 0, 0, 0];       % normalized Gell-Mann matrix lambda 1
X2 = [0, 0, c_GM; 0, 0, 0; c_GM, 0, 0];      % lambda 4
X3 = [0, 0, 0; 0, 0, c_GM; 0, c_GM, 0];      % lambda 6
Y1 = [0, -1i*c_GM, 0; 1i*c_GM, 0, 0; 0, 0, 0];   % lambda 2
Y2 = [0, 0, -1i*c_GM; 0, 0, 0; 1i*c_GM, 0, 0];   % lambda 5
Y3 = [0, 0, 0; 0, 0, -1i*c_GM; 0, 1i*c_GM, 0];   % lambda 7
Z1 = [c_GM, 0, 0; 0, -c_GM, 0; 0, 0, 0];      % lambda 3
Z2 = [c_GM/(3)^.5, 0, 0; 0, c_GM/(3)^.5, 0; 0, 0, -2*c_GM/(3)^.5];      % lambda 8

M1 = X1; % x axis
M2 = Z2; % y axis

A = zeros(dA,dA,N,k);   % children elements

% Generate childrens
for i=1:N
    A(:,:,i,1) = (1/dA)*(eye(dA)) + (1/2)*(cos((i-1)*(2*pi)/(2*N)) * M1 + sin((i-1)*(2*pi)/(2*N)) * M2)    %#ok<NOPTS> 
    A(:,:,i,2) = (1/dA)*eye(dA)   %#ok<NOPTS> 
    A(:,:,i,3) = eye(dA) - A(:,:,i,1) - A(:,:,i,2)    %#ok<NOPTS> 
end 

[F_ax, eta, flag] = WNR_sdp1(A,rho_AB);   % 10(a) in paper

eta    %#ok<NOPTS> 

