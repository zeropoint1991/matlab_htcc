%# cat speedup.m

function s = speedup(n)

a = rand(n);
b = rand(n);
tic;
c = a*b;
time1 = toc;

ga = gpuArray(a);
gb = gpuArray(b);
tic;
gc = ga*gb;
time2 = toc;

s = time1/time2;