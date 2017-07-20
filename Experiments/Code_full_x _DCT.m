clc; clear all; close all;
tic;
rng('default');

Img_name = 'einstein.pgm'
% Img_name = 'lena.pgm'
% Img_name = 'man.tiff'
% Img_name = 'moon.tif'
I =  double(imread(Img_name)); % reading the image
I = imresize(I,[128 128]); % resizing the image

[m,n] = size(I);

%qmf = MakeONFilter('Coiflet',2);
%X = FWT2_PO(I,3,qmf);
X = dct2(I); % converting the image data to the sparse domain
x = X(:); % vectorizing the sparse data

N = n*n; % length of the signal
M = 4000; % No. of measurements
%q = randperm(N)';
%g = randperm(N,M);

for K = 1000:200:2000
    K
    
    A = randn(M,N);
    for i = 1:N % normalizing the coloumns of measurement matrix
        r = A(:,i);
        nrm1 = norm(r);
        A(:,i) = r/nrm1;
    end
    
    b = A*x;
    
    [Xest,support] = SP(A,b,K); % calling SP function
    
    
    rs1 = reshape(Xest,[128 128]);
    inv_dct = idct2(rs1);
    %inv_dct = IWT2_PO(rs1,3,qmf);
    mseImg = (I-double(inv_dct)).^2;
    MSE = sum(mseImg(:))/(m*n)
    PSNR = 10*log10(255^2/MSE)
end
toc;

