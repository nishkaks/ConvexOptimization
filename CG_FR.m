function [inform,x] = CG_FR(fun,x,sdparams)

k = 1;
alpha = 1;
params = struct('c1',0.01,'c2',0.3,'maxit',100);
x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
pk = -x.g;

while (norm(x.g,inf) > (sdparams.toler * (1+abs(x.f))) && k <= sdparams.maxit)   
    gradfk = x.g;    
    [alpha,x] = StepSize(fun, x, pk, alpha, params);
    betakplus1 = (x.g' * x.g)/(gradfk' * gradfk);
    pk = -x.g + (betakplus1 * pk);    
    k= k+1;
end;

inform.iter = k-1;
if (norm(x.g,inf) <= (sdparams.toler * (1+abs(x.f))))   
    inform.status = 1;
else
    inform.status = 0;
end
