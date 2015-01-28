function [alfa,x] = StepSize(fun, x, d, alfa, params)
alphaMax = realmax;

i = 1;
x0 = x.p;
alphaiMinus1 = 0;
alphai = alfa;
phi.Zero = x.f;
gradPhi.Zero = x.g' * d;
phi.AlphaiMinus1 = phi.Zero;
gradPhi.AlphaiMinus1 = gradPhi.Zero;

while (alphai < alphaMax && i < params.maxit)
    x.p = x0 + alphai * d;
    x.f = feval(fun,x.p,1);    
    phi.Alphai = x.f;
    if((phi.Alphai > phi.Zero + (params.c1 * alphai * gradPhi.Zero)) || ( (i > 1) && phi.Alphai >= phi.AlphaiMinus1))
       %alphaLo = alphaiMinus1;
       %alphaHi = alphai;
       [alfa,x] = Zoom(fun,x0,d,alphaiMinus1,alphai,params,phi.AlphaiMinus1,gradPhi.AlphaiMinus1,phi.Alphai,phi.Zero,gradPhi.Zero);
       return;
    end
    
    x.g = feval(fun,x.p,2);
    gradPhi.Alphai = x.g'*d;
    if (abs(gradPhi.Alphai) <= params.c2 * abs(gradPhi.Zero) )
        alfa = alphai;
        return;
    end
    
    if (gradPhi.Alphai >= 0)
       %alphaHi = alphaiMinus1;
       %alphaLo = alphai;    
       [alfa,x] = Zoom(fun,x0,d,alphai,alphaiMinus1,params,phi.Alphai,gradPhi.Alphai,phi.AlphaiMinus1,phi.Zero,gradPhi.Zero);
       return;
    end  
    
    alphaiMinus1 = alphai;
    alphai = 3 * alphaiMinus1;
    phi.AlphaiMinus1 = phi.Alphai;
    gradPhi.AlphaiMinus1 = gradPhi.Alphai;
    i= i+1;    
end;


function [alpha,x] = Zoom(fun,x0,d,alphaLo,alphaHi,params,phiLo,gradPhiLo,phiHi,phiZero,gradPhiZero)
j=1;
while (j < params.maxit)
    alphaDiff = alphaHi - alphaLo;
    alphaJ = alphaLo - ((gradPhiLo * alphaDiff^2)/(2* (phiHi - phiLo - gradPhiLo * alphaDiff)));
    safeguardLo = min(alphaLo,alphaHi) + 0.2 * abs(alphaDiff);
    safeguardHi = max(alphaLo,alphaHi) - 0.2 * abs(alphaDiff);
    if (alphaJ < safeguardLo)
        alphaJ = safeguardLo; 
    elseif (alphaJ >  safeguardHi)
        alphaJ = safeguardHi;    
    end
    xj = x0 + alphaJ * d;
    phiAlphaj = feval(fun,xj,1);
    if((phiAlphaj > phiZero + (params.c1 * alphaJ * gradPhiZero)) || (phiAlphaj >= phiLo))
        alphaHi = alphaJ;
        phiHi = phiAlphaj;
    else
        xg = feval(fun,xj,2);
        gradPhiAlphaj = xg'  * d;
        if (abs(gradPhiAlphaj) <= params.c2 * abs(gradPhiZero) )
            alpha = alphaJ; 
            x.p = xj;
            x.f = phiAlphaj;
            x.g = xg;
            return;
        end
        if (gradPhiAlphaj * (alphaDiff) >= 0)
            alphaHi = alphaLo;
            phiHi = phiLo;
        end
        alphaLo = alphaJ;
        phiLo = phiAlphaj;
        gradPhiLo = gradPhiAlphaj;
    end
    j= j+1;
end



