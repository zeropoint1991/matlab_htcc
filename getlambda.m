function lambda = getlambda(popsize)
%2Œ¨œÚ¡ø
lambda=ones(popsize,2);
lambda(:,1)=(1:popsize)/popsize;
lambda(:,2)=lambda(:,2)-lambda(:,1);
end

