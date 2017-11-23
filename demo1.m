%# cat demo1.m

for i=1000:2000:9000
disp(sprintf('GPU speedup for array multiplicaton of
array size %d is %.0f', i, speedup(i)));
end