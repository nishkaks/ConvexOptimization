function [inform,x] = SteepDescent(fun,x,sdparams)

k = 1;
alpha = 1;
params = struct('c1',0.01,'c2',0.3,'maxit',100);
x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);

while (norm(x.g,inf) > (sdparams.toler * (1+abs(x.f))) && k <= sdparams.maxit)    
    [alpha,x] = StepSize(fun, x, -1 * x.g, 1.2 * alpha, params);
    k= k+1;
end;

inform.iter = k-1;
if (norm(x.g,inf) <= (sdparams.toler * (1+abs(x.f))))   
    inform.status = 1;
else
    inform.status = 0;
end
