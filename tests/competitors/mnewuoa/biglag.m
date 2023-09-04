function [D,alpha] = biglag(n,npt,Xopt,Xpt,Bmat,Zmat,idz,knew,delta)

% n is the number of variables.
% npt is the number of interpolation equations.
% Xopt is the best interpolation point so far.
% Xpt contains the coordinates of the current interpolation points.
% Bmat provides the last n columns of H.
% Zmat and idz give a factorization of the first npt by npt submatrix of H.
% knew is the index of the interpolation point that is going to be moved.
% delta is the current trust region bound.
% D will be set to the step from Xopt to the new point.
% alpha will be set to the knew-th diagonal element of the H matrix.
%
% The step D is calculated in a way that attempts to maximize the modulus
% of LFUNC(Xopt+D), subject to the bound ||D|| .LE. delta, where LFUNC is
% the knew-th Lagrange function.

% Set some constants.
delsq = delta*delta;

% Set the first npt components of Hcol to the leading elements of the
% knew-th column of H.
iterc = 0;

%------------------
Hcol = -Zmat(knew,1:idz-1)*Zmat(:,1:idz-1)';
    Hcol = Hcol +Zmat(knew,idz:npt-n-1)*Zmat(:,idz:npt-n-1)';
% Hcol = zeros(1,npt);
% for j=1:npt-n-1
%     temp = Zmat(knew,j);
%     if (j<idz)
%         temp = -temp;
%     end
%     Hcol = Hcol +temp*Zmat(:,j)';
% end
%-------------------

alpha = Hcol(knew);

% Set the unscaled initial direction D. Form the gradient of LFUNC at
% Xopt, and multiply D by the second derivative matrix of LFUNC.
D = Xpt(knew,:)-Xopt;
dd = sum(D.^2);

%------------------
Gc = Bmat(knew,:)+(Hcol.*(Xopt*Xpt'))*Xpt;
Gd = (Hcol.*(D*Xpt'))*Xpt;
% Gc = Bmat(knew,:);
% Gd = zeros(1,n);
% for k=1:npt
%     temp = Hcol(k)*(Xpt(k,:)*Xopt');
%     tsum = Hcol(k)*(Xpt(k,:)*D');
%     Gc = Gc+temp*Xpt(k,:);
%     Gd = Gd+tsum*Xpt(k,:);
% end
%------------------

% Scale D and Gd, with a sign change if required. Set S to another
% vector in the initial two dimensional subspace.
gg = sum(Gc.^2);
sp = D*Gc';
dhd = D*Gd';
scale = delta/sqrt(dd);
if (sp*dhd<0)
    scale = -scale;
end
temp = 0;
if (sp*sp>.99*dd*gg)
    temp = 1;
end
tau = scale*(abs(sp)+.5*scale*abs(dhd));
if (gg*delsq<.01*tau*tau)
    temp = 1;
end
D = scale*D;
Gd = scale*Gd;
S = Gc + temp*Gd;

while (iterc<n)
    % Begin the iteration by overwriting S with a vector that has the
    % required length and direction, except that termination occurs if
    % the given D and S are nearly parallel.
    iterc = iterc+1;
    dd = sum(D.^2);
    sp = D*S';
    ss = sum(S.^2);
    temp = dd*ss-sp*sp;
    if (temp<=(1e-8)*dd*ss)
        return;
    end
    denom = sqrt(temp);
    S = (dd*S-sp*D)/denom;

    % Calculate the coefficients of the objective function on the circle,
    % beginning with the multiplication of S by the second derivative matrix.
    
    %------------------
    W = (Hcol.*(S*Xpt')) * Xpt;
%     W = zeros(1,n);
%     for k=1:npt
%         tsum = Hcol(k)*(Xpt(k,:)*S');
%         W = W+tsum*Xpt(k,:);
%     end
    %------------------
    cf1 = .5*(S*W');
    cf2 = D*Gc';
    cf3 = S*Gc';
    cf4 = .5*(D*Gd')-cf1;
    cf5 = S*Gd';

    % Seek the value of the angle that maximizes the modulus of tau.
    taubeg = cf1+cf2+cf4;
    taumax = taubeg;
    tauold = taubeg;
    isave = 0;
    iu = 49;
    temp = 2*pi/(iu+1);
    for i=1:iu
        angle = i*temp;
        cth = cos(angle);
        sth = sin(angle);
        tau = cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth;
        if (abs(tau)>abs(taumax))
            taumax = tau;
            isave = i;
            tempa = tauold;
        elseif (i==isave+1)
            tempb = tau;
        end
        tauold = tau;
    end
    if (isave==0)
        tempa = tau;
    end
    if (isave==iu)
        tempb = taubeg;
    end
    step = 0;
    if (tempa~=tempb)
        tempa = tempa-taumax;
        tempb = tempb-taumax;
        step = .5*(tempa-tempb)/(tempa+tempb);
    end
    angle = temp*(isave+step);

    % Calculate the new D and Gd. Then test for convergence.
    cth = cos(angle);
    sth = sin(angle);
    tau = cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth;
    D = cth*D+sth*S;
    Gd = cth*Gd+sth*W;
    S = Gc+Gd;

    if (abs(tau)<=1.1*abs(taubeg))
        return;
    end
end
return;
