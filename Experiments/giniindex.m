function[gi] = giniindex(x)

[N,M] = size(x);
fk = 0;


f = sort(abs(x),'ascend');

for i = 1:N
    
    fk = fk + (f(i)/norm(f,1)) * ((N - i + 0.5)/N);
    %fk = fk + (f(i)/(sqrt(N)*norm(x))) * ((N - i + 0.5)/N);
    
end

gi = 1 - 2*fk;
