function u = Anti_TV(g,mu,gamma)
%function u = SB_ATV(g,mu,gamma)
% Split Bregman Anisotropic Total Variation Denoising
% This function performs TV norm based optimization for contrast 
% enhancementon the input image using the given 'mu' value and 'gamma' 
% value as given in the paper

%Ravi, H., Subramanyam, A. V., Emmanuel, S., "ACE - An Effective
%Anti-forensic Contrast Enhancement Technique", accepted in IEEE Signal
%Processing Letters (IEEE SPL), 2015. Please cite the paper if you use this
%code.

%Input
% g - input image
% mu - value used as given in eq 4 in the paper
% gamma - value for contrast enhancement
%Output
% u - Output image that is enhanced anti-forensically

%This code might not be optimized and completely free of errors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Code written by, @ Hareesh Ravi (Research Associate at IIITD)    %%%%  
%%%%  (haree.24@gmail.com)                                             %%%% 
%%%%  code can be used and modified for research purposes.             %%%% 
%%%%  Kindly let me know of mistakes by mail                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initialize normal enhanced image as in second term of the optimization
%equation
gHS = uint8(imadjust(g));
gGC = uint8(255.*((double(g)./255).^(gamma)));

%Initialize other stuff based on the terms used in eq 4 in the paper
g = double(g(:));
n = length(g);
b = zeros(2*n,1);
d = b;
u = g;
err = 1;k = 1;
tol = 1e-3;
lambda = 0.05;

%Initialize the d_x and d_y used in the equation 4 in the paper. For now it
%is only for square input images, but can be modified. 
[B, Bt, BtB] = DiffOper(sqrt(n));

%Run iterations until tolerance is reached
while err > tol
    up = u;
    %Perforn conjugated gradient descent using inbuilt function
    if gamma == 100
        [u,~] = cgs((1.6.*speye(n)+lambda.*BtB), ((0.2.*g)+(1.4.*double(gHS(:)))-lambda*Bt*(b-d)),1e-5,100);
    elseif gamma>1        
        [u,~] = cgs((1.6.*speye(n)+lambda.*BtB), (0.4.*g+1.2.*double(gGC(:))-lambda*Bt*(b-d)),1e-5,100); 
    else
        [u,~] = cgs((1.6.*speye(n)+lambda.*BtB), (0.1.*g+1.5.*double(gGC(:))-lambda*Bt*(b-d)),1e-5,100); 
    end
    
    %Perfom shrinkage operations as given in the paper
    Bub = B*u+b;
    d = max(abs(Bub)-mu/lambda,0).*sign(Bub);
    b = Bub-d;
    
    %Calculate error
    err = norm(up-u)/norm(u);
    %fprintf('err=%g \n',err);
    k = k+1;
end
%fprintf('Stopped because norm(up-u)/norm(u) <= tol=%.1e\n',tol);
end

function [B, Bt, BtB] = DiffOper(N)
%function [B, Bt, BtB] = DiffOper(N)
%This function generates the difference operators used in the TV norm
%equation as described in the SPL paper. Check paper for more details.
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;
end
