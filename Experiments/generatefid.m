clc; clear all; close all;

load('Biomed_DWT_baomp.mat')

for i = 1:2
    k = StoreforDiffM(:,i);
    rs1 = reshape(k,[64 64]);
    inv_dct = IWT2_PO(rs1,3,qmf);
    %inv_dct = idct2(rs1);
    chngSize = imresize(inv_dct,[256 256]);
    figure,imshow(chngSize,[])
    
end