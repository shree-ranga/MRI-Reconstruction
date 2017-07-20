
% baomp + giniindex final

% Author: Shree Ranga Raju NM
% SSP Lab
% Dept of ECE, Indian Institute of Science


function [xest, estsupport] = bagiomp(A,b,maxiter)

[m,n] = size(A); %size of the measurement matrix

finalSet = []; % finalist set

deleteSetor = []; % set to delete wrongly chosen atoms

candidateSet = [];





r = b; % initial residue

maxerror = 10^-6; % maximum error

for i = 1:maxiter
    
    
    trans = abs(A'*r); % matched filter
    
    trans(finalSet) = 0; % not interested in these locations
    
    maxval1 = max(trans); % maximum value of the matched filter
    
    ginir = giniindex(r);
    
    threshold1 = ginir*maxval1; % caluclating threshold using Gini Index
    
    
    index = find(trans >= threshold1); % finding all the indices which are greater than the threshold value
    
    candidateSet = index'; % forming the candidate set
    
    unionset = [candidateSet,finalSet]; % union of the candidate set and the final set
    
    xest = A(:,unionset)\b; % solving the least square problem
    
    maxval2 = max(abs(xest(1:length(candidateSet)))); % finding the maximum absolute value of xest
    
    giXest = giniindex(xest);
    
    threshold2 = giXest*maxval2 ;% caluclating threshold using Gini Index
    
    
    
    threshold2Elements = find(abs(xest) < threshold2); % finding all the elements of xest which falls below the threshold
    
    deleteSetor = unionset(threshold2Elements); % selects the indices corresponding to the threshold elements from the union set
    
    finalSet = setdiff(unionset,deleteSetor); % formation of the final support set
    
    xest = A(:,finalSet)\b; % caluclating the least square problem again
    
    r = b - A(:,finalSet)*xest; % calulclating the residue
    
    
    nrmr = norm(r);
    
    
    if(nrmr < maxerror)
        break;
    end
    
end

estsupport = finalSet; % final support set of the signal

xest = zeros(n,1);

xest(estsupport) = A(:,finalSet)\b; % approximate of the true signal



