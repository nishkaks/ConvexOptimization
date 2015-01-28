global numf numg 

% declare dimension
n = 100;
n4 = n/4;

% set up initial point
z = [ 3 -1 0 1]'; ztemp = ones(n4,1);
x0 = kron(ztemp,z);
x=struct('p',x0);

% solve with Fletcher-Reeves CG

nonCGparams = struct('maxit',10000,'toler',1.0e-5)  ;
fprintf(' \n** Fletcher-Reeves CG on xpowsing\n');
numf=0; numg=0;
[inform,x] = CG_FR(@xpowsing,x,nonCGparams);
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %8.4g.\n', mNparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('\n  Ending value: '); fprintf('%8.4g ',x.f);
fprintf('; No. function evaluations: %d',numf);
fprintf('; No. gradient evaluations %d',numg);
fprintf('\n  Norm of ending gradient: %8.4g\n\n', norm(x.g,inf));


% solve with Polak-Ribiere CG with modification

x=struct('p',x0);
nonCGparams = struct('maxit',10000,'toler',1.0e-5)  ;
fprintf(' \n** PR+ CG on xpowsing\n');
numf=0; numg=0;
[inform,x] = CG_PRplus(@xpowsing,x,nonCGparams);
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %8.4g.\n', mNparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('\n  Ending value: '); fprintf('%8.4g ',x.f);
fprintf('; No. function evaluations: %d',numf);
fprintf('; No. gradient evaluations %d',numg);
fprintf('\n  Norm of ending gradient: %8.4g\n\n', norm(x.g,inf));


% solve with steepest descent, with the same line search

x=struct('p',x0);
sdparams = struct('maxit',10000,'toler',1.0e-5);
fprintf(' \n** Steepest Descent on xpowsing\n');
numf=0; numg=0;
[inform,x] = SteepDescent(@xpowsing,x,sdparams);
if inform.status == 0
  fprintf('CONVERGENCE FAILURE: %d steps were taken without\n', inform.iter);
  fprintf('gradient size decreasing below %8.4g.\n', sdparams.toler);
else
  fprintf('Success: %d steps taken\n', inform.iter);
end
fprintf('\n  Ending value: '); fprintf('%8.4g ',x.f);
fprintf('; No. function evaluations: %d',numf);
fprintf('; No. gradient evaluations %d',numg);
fprintf('\n  Norm of ending gradient: %8.4g\n\n', norm(x.g,inf));
