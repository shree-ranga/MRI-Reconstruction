% Gini Index based BAOMP algorithm
% Author : Shree Ranga Raju N M
% SYDE 633 - Remote Sensing Systems Project Code

clear all; close all; clc
%rng('default');
%
Img = imread('mrImage.jpg');
Img = rgb2gray(Img);
I = double(imresize(Img,[64 64]));
%figure, imshow(I,[])
[m n] = size(I);

N = n*n;
%qmf = MakeONFilter('Coiflet',2);
%X = FWT2_PO(I,3,qmf);
X = dct2(I);
x = X(:);


M = round(N*[0.3,0.6]);

for j = 1: length(M)
    
    A = randn(M(j),N);
    for i = 1:N % normalizing the coloumns of measurement matrix
        r = A(:,i);
        nrm1 = norm(r);
        A(:,i) = r/nrm1;
    end
    disp('Normalizing of measurement matrix done!')
    
    b = A*x;
    
    
    [xest,support] = bagiomp(A,b,M(j));
    
    disp('Calculated the estimated signal and support set');
    
    StoreforDiffM(:,j) = xest;
    rs1 = reshape(xest,[64 64]);
    %inv_dct = IWT2_PO(rs1,3,qmf);
    inv_dct = idct2(rs1);
    figure,imshow(inv_dct,[])
    mseImg = (I-double(inv_dct)).^2;
    MSE = sum(mseImg(:))/(m*n);
    PSNR(j) = 10*log10(255^2/MSE)
    
    
end

save('Biomed_DCT_bagiomp.mat')


