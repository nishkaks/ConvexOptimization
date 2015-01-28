mu=0.01; L=1; kappa=L/mu;
n=100;
A = randn(n,n); [Q,R]=qr(A);
D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D = mu + D*(L-mu);
A = Q'*diag(D)*Q;
epsilon=1.e-6;
kmax=2000;

hold on;

klist = zeros(1,10);
%averagerate = zeros(1,10);

% use the same set of random starts for all the algorithms
x0list = zeros(n,10);
for i = 1:1:10
    x0list(:,i) = randn(n,1);
end

%--------------------------------------%
%  Steepest Descent - fixed steps      %
%--------------------------------------%
alpha = 2/(mu + L);
for i = 1:1:10    
    x0 = x0list(:,i); 
    
    k = 1;
    fk = 0.5 * x0' * A * x0;
    fplot = fk;
    xk = x0;    
    %rate = 0;    
    while (fk >= epsilon && k < kmax)
        k= k+1;
        grad = (A) * xk;
        xk1 = xk - alpha * grad; 
        fk = 0.5 * xk1' * A * xk1;
        xk = xk1; 
        fplot = [fplot;fk];
        %newrate = (fplot(k) / fplot(k-1));
        %if newrate > rate
        %    rate = newrate;
        %end;
    end;
    
    %averagerate(i) = rate;
    
    klist(i) = k;
    
    y = log(fplot);
    x = 1:1:k;    
    plot(x,y,'r')
    %hold on;
    
end
av_sd = sum(klist) / 10;
fprintf(1,' steepest descent - fixed steps : %7.1f\n', av_sd);
%avgrate_sd = mean(averagerate);
%fprintf(1,' steepest descent - fixed steps, rate: %7.1f\n', avgrate_sd);

%--------------------------------------%
% Steepest Descent - exact line search %
%--------------------------------------%
klist = zeros(1,10);
%averagerate = zeros(1,10);
for i = 1:1:10    
    x0 = x0list(:,i);
    
    k = 1;
    fk = 0.5 * x0' * A * x0;
    fplot = fk;
    xk = x0;    
    %rate = 0;    
    while (fk >= epsilon && k < kmax)
        k= k+1;
        grad = (A) * xk;
        alpha = (xk' * A * A * xk)/(xk' * A * A * A * xk);
        xk1 = xk - alpha * grad; 
        fk = 0.5 * xk1' * A * xk1;
        xk = xk1; 
        fplot = [fplot;fk];
        %newrate = (fplot(k) / fplot(k-1));
        %if newrate > rate
        %    rate = newrate;
        %end;
    end;
    
    %averagerate(i) = rate;
    
    klist(i) = k;
    
    y = log(fplot);
    x = 1:1:k;    
    plot(x,y,'g')
    %hold on;
    
end
av_sde = sum(klist) / 10;
fprintf(1,' steepest descent - exact steps : %7.1f\n', av_sde);
%avgrate_sde = mean(averagerate);
%fprintf(1,' steepest descent - exact steps rate: %7.1f\n', avgrate_sde);


%--------------------------------------%
%           Heavy Ball                 %
%--------------------------------------%
alpha = (4/L) * (1/(1+(1/sqrt(kappa)))^2);
beta = (1 - (2/(sqrt(kappa) + 1)))^2;
klist = zeros(1,10);
%averagerate = zeros(1,10);
for i = 1:1:10    
    x0 = x0list(:,i);
    
    k = 1;
    fk = 0.5 * x0' * A * x0;
    fplot = fk;
    xk = x0;    
    %rate = 0;  
    
    % Need two iterations for heavy ball
    k = k+1;
    grad = (A) * xk;
    xk1 = xk - alpha * grad; 
    fk = 0.5 * xk1' * A * xk1;
    fplot = [fplot;fk];
    xkminus1 = xk;
    xk = xk1;
    
    while (fk >= epsilon && k < kmax)
        k= k+1;
        grad = (A) * xk;
        xk1 = xk - alpha * grad + beta * (xk - xkminus1); 
        fk = 0.5 * xk1' * A * xk1;
        xkminus1 = xk;
        xk = xk1; 
        fplot = [fplot;fk];
        %newrate = (fplot(k) / fplot(k-1));
        %if newrate > rate
        %    rate = newrate;
        %end;
    end;
    
    %averagerate(i) = rate;
    
    klist(i) = k;
    
    y = log(fplot);
    x = 1:1:k;    
    plot(x,y,'b')
    
    
end
av_hb = sum(klist) / 10;
fprintf(1,' heavy ball : %7.1f\n', av_hb);
%avgrate_hb = mean(averagerate);
%fprintf(1,' heavy ball rate: %7.1f\n', avgrate_hb);

%--------------------------------------%
%           Conjugate gradient         %
%--------------------------------------%

klist = zeros(1,10);
%averagerate = zeros(1,10);
for i = 1:1:10    
    x0 = x0list(:,i);    
    k = 1;
    
    rk = (A * x0);
    pk = - rk;
    fk = 0.5 * x0' * A * x0;
    fplot = fk;
    xk = x0;    
    %rate = 0; 
    while (fk >= epsilon && k < kmax) % not used any(rk >= epsilon) 
                                      % for comparison
        
        alpha = -1 * (rk' * pk)/(pk' * A * pk);
        xk1 = xk + alpha * pk; 
        fk = 0.5 * xk1' * A * xk1;
        fplot = [fplot;fk];
        xk = xk1;
        
        rk = A * xk1;
        beta = (rk' * A * pk)/(pk' * A * pk);
        pk1 = - rk + (beta * pk);
        pk = pk1;
        
        k= k+1;
       
        %newrate = (fplot(k) / fplot(k-1));
        %if newrate > rate
        %    rate = newrate;
        %end;
    end;
   
    %averagerate(i) = rate;
    
    klist(i) = k;
    
    y = log(fplot);
    x = 1:1:k;    
    plot(x,y,'m')
    %hold on;    
end

av_cg = sum(klist) / 10;
fprintf(1,' conjugate gradient : %7.1f\n', av_cg);
%avgrate_cg = mean(averagerate);
%fprintf(1,' conjugate gradient,rate : %7.1f\n', avgrate_cg);

%---------------------------------------------------%
% Hybrid - steepest descent - fixed and Heavy Ball  %
%---------------------------------------------------%

beta = (1 - (2/(sqrt(kappa) + 1)))^2;
klist = zeros(1,10);
%averagerate = zeros(1,10);
for i = 1:1:10    
    x0 = x0list(:,i);
    
    k = 1;
    fk = 0.5 * x0' * A * x0;
    fplot = fk;
    xk = x0;    
    %rate = 0;  
    
    % Need two iterations for heavy ball
    k = k+1;
    grad = (A) * xk;
    xk1 = xk - alpha * grad; 
    fk = 0.5 * xk1' * A * xk1;
    fplot = [fplot;fk];
    xkminus1 = xk;
    xk = xk1;
    alpha = 2/(mu + L);
    
    while (fk >= epsilon && k < kmax)
        k= k+1;
        grad = (A) * xk;
        if (k>25)
           alpha = (4/L) * (1/(1+(1/sqrt(kappa)))^2);
        end;
        xk1 = xk - alpha * grad + beta * (xk - xkminus1); 
        fk = 0.5 * xk1' * A * xk1;
        xkminus1 = xk;
        xk = xk1; 
        fplot = [fplot;fk];
        %newrate = (fplot(k) / fplot(k-1));
        %if newrate > rate
        %    rate = newrate;
        %end;
    end;
    %averagerate(i) = rate;
    
    klist(i) = k;
    
    y = log(fplot);
    x = 1:1:k;    
    plot(x,y,'c')
    
    
end
av_sdhb = sum(klist) / 10;
fprintf(1,' steepest descent - fixed + heavy ball : %7.1f\n', av_sdhb);
%avgrate_sdhb = mean(averagerate);
%fprintf(1,' steepest descent - fixed + heavy ball rate: %7.1f\n', avgrate_sdhb);
