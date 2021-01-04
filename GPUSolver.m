function [x] = GPUSolver(A,b)
x=gpuArray(gather(A) \ gather(b));
end

