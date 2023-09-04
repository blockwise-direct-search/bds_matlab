function [D,W,Vlag,beta] = bigden(n,npt,Xopt,Xpt,Bmat,Zmat,idz,kopt,knew,D)
% n is the number of variables.
% npt is the number of interpolation equations.
% Xopt is the best interpolation point so far.
% Xpt contains the coordinates of the current interpolation points.
% Bmat provides the last n columns of H.
% Zmat and idz give a factorization of the first npt by npt submatrix of H.
% kopt is the index of the optimal interpolation point.
% knew is the index of the interpolation point that is going to be moved.
% D will be set to the step from Xopt to the new point, and on entry it
%   should be the D that was calculated by the last call of BIGLAG. The
%   length of the initial D provides a trust region bound on the final D.
% W will be set to Wcheck for the final choice of D.
% Vlag will be set to Theta*Wcheck+e_b for the final choice of D.
% beta will be set to the value that will occur in the updating formula
%   when the knew-th interpolation point is moved to its new position.
%
% D is calculated in a way that should provide a denominator with a large
% modulus in the updating formula when the knew-th interpolation point is
% shifted to the new position Xopt+D.

% Store the first npt elements of the knew-th column of H in W(n+1)
% to W(n+npt).
W(n+[1:npt]) = -Zmat(knew,1:idz-1)*Zmat(:,1:idz-1)';
W(n+[1:npt]) = W(n+[1:npt])+Zmat(knew,idz:npt-n-1)*Zmat(:,idz:npt-n-1)';
alpha = W(n+knew);

% The initial search direction D is taken from the last call of BIGLAG,
% and the initial S is set below, usually to the direction from X_OPT
% to X_KNEW, but a different direction to an interpolation point may
% be chosen, in order to prevent S from being nearly parallel to D.
dd = sum(D.^2);
S = Xpt(knew,1:n)-Xopt;
ds = D*S';
ss = sum(S.^2);
xoptsq = sum(Xopt.^2);
if (ds*ds>.99*dd*ss)
    ksav = knew;
    dtest = ds*ds/ss;
    for k=1:npt
        if (k~=kopt)
            dstemp = D*(Xpt(k,1:n)-Xopt)';
            sstemp = sum((Xpt(k,1:n)-Xopt).^2);
            if (dstemp*dstemp/sstemp<dtest)
                ksav = k;
                dtest = dstemp*dstemp/sstemp;
                ds = dstemp;
                ss = sstemp;
            end
        end
    end
    S = Xpt(ksav,1:n)-Xopt;
end
ssden = dd*ss-ds*ds;
iterc = 0;
densav = 0;

% Begin the iteration by overwriting S with a vector that has the
% required length and direction.
while (1)
    iterc = iterc+1;
    temp = 1/sqrt(ssden);
    S = temp*(dd*S-ds*D);
    xoptd = Xopt*D';
    xopts = Xopt*S';

    % Set the coefficients of the first two terms of beta.
    tempa = .5*xoptd*xoptd;
    tempb = .5*xopts*xopts;
    Den = zeros(1,9);
    Den(1) = dd*(xoptsq+.5*dd)+tempa+tempb;
    Den(2) = 2*xoptd*dd;
    Den(3) = 2*xopts*dd;
    Den(4) = tempa-tempb;
    Den(5) = xoptd*xopts;

    % Put the coefficients of Wcheck in Wvec.
    Wvec = zeros(npt+n,5);
    %------------------
    %     for k=1:npt
    %         tempa = Xpt(k,1:n)*D';
    %         tempb = Xpt(k,1:n)*S';
    %         tempc = Xpt(k,1:n)*Xopt';
    %         Wvec(k,1) = .25*(tempa*tempa+tempb*tempb);
    %         Wvec(k,2) = tempa*tempc;
    %         Wvec(k,3) = tempb*tempc;
    %         Wvec(k,4) = .25*(tempa*tempa-tempb*tempb);
    %         Wvec(k,5) = .5*tempa*tempb;
    %     end
    tempa = Xpt(:,1:n)*D';
    tempb = Xpt(:,1:n)*S';
    tempc = Xpt(:,1:n)*Xopt';
    Wvec(1:npt,1) = .25*(tempa.^2+tempb.^2);
    Wvec(1:npt,2) = tempa.*tempc;
    Wvec(1:npt,3) = tempb.*tempc;
    Wvec(1:npt,4) = .25*(tempa.^2-tempb.^2);
    Wvec(1:npt,5) = .5*tempa.*tempb;
    %------------------
    
    Wvec(npt+[1:n],2) = D;
    Wvec(npt+[1:n],3) = S;

    % Put the coefficents of THETA*Wcheck in Prod.
    Prod = zeros(npt+n,5);
    for jc=1:5
        nw = npt;
        if (jc==2) || (jc==3)
            nw = npt+n;
        end

        %------------------
        Prod(1:npt,jc) = -Zmat(:,1:idz-1)*(Zmat(:,1:idz-1)'*Wvec(1:npt,jc));
        Prod(1:npt,jc) = Prod(1:npt,jc)+Zmat(:,idz:npt-n-1)*(Zmat(:,idz:npt-n-1)'*Wvec(1:npt,jc));
%         for j=1:npt-n-1
%             tsum = Zmat(1:npt,j)'*Wvec(1:npt,jc);
%             if (j<idz)
%                 tsum = -tsum;
%             end
%             Prod(1:npt,jc) = Prod(1:npt,jc)+tsum*Zmat(1:npt,j);
%         end
        %--------------------

        if (nw==npt+n)
            Prod(1:npt,jc) = Prod(1:npt,jc)+Bmat(1:npt,1:n)*Wvec(npt+[1:n],jc);
        end
        %--------------------
        Prod(npt+[1:n],jc) = Bmat(1:nw,1:n)'*Wvec(1:nw,jc);
%             
%         for j=1:n
%             tsum = Bmat(1:nw,j)'*Wvec(1:nw,jc);
%             Prod(npt+j,jc) = tsum;
%         end
        %norm(Prod1(npt+[1:n],jc)-Prod(npt+[1:n],jc))
    end

    % Include in Den the part of beta that depends on THETA.
    for k=1:npt+n
        Par(1:5) = .5*Wvec(k,1:5)*Prod(k,1:5)';
        tsum = sum(Par);
        Den(1) = Den(1)-Par(1)-tsum;
        tempa = Prod(k,1)*Wvec(k,2)+Prod(k,2)*Wvec(k,1);
        tempb = Prod(k,2)*Wvec(k,4)+Prod(k,4)*Wvec(k,2);
        tempc = Prod(k,3)*Wvec(k,5)+Prod(k,5)*Wvec(k,3);
        Den(2) = Den(2)-tempa-.5*(tempb+tempc);
        Den(6) = Den(6)-.5*(tempb-tempc);
        tempa = Prod(k,1)*Wvec(k,3)+Prod(k,3)*Wvec(k,1);
        tempb = Prod(k,2)*Wvec(k,5)+Prod(k,5)*Wvec(k,2);
        tempc = Prod(k,3)*Wvec(k,4)+Prod(k,4)*Wvec(k,3);
        Den(3) = Den(3)-tempa-.5*(tempb-tempc);
        Den(7) = Den(7)-.5*(tempb+tempc);
        tempa = Prod(k,1)*Wvec(k,4)+Prod(k,4)*Wvec(k,1);
        Den(4) = Den(4)-tempa-Par(2)+Par(3);
        tempa = Prod(k,1)*Wvec(k,5)+Prod(k,5)*Wvec(k,1);
        tempb = Prod(k,2)*Wvec(k,3)+Prod(k,3)*Wvec(k,2);
        Den(5) = Den(5)-tempa-.5*tempb;
        Den(8) = Den(8)-Par(4)+Par(5);
        tempa = Prod(k,4)*Wvec(k,5)+Prod(k,5)*Wvec(k,4);
        Den(9) = Den(9)-.5*tempa;
    end

    % Extend Den so that it holds all the coefficients of DENOM.
    Par(1:5) = .5*Prod(knew,1:5).^2;
    tsum = sum(Par);
    Denex(1) = alpha*Den(1)+Par(1)+tsum;
    tempa = 2*Prod(knew,1)*Prod(knew,2);
    tempb = Prod(knew,2)*Prod(knew,4);
    tempc = Prod(knew,3)*Prod(knew,5);
    Denex(2) = alpha*Den(2)+tempa+tempb+tempc;
    Denex(6) = alpha*Den(6)+tempb-tempc;
    tempa = 2*Prod(knew,1)*Prod(knew,3);
    tempb = Prod(knew,2)*Prod(knew,5);
    tempc = Prod(knew,3)*Prod(knew,4);
    Denex(3) = alpha*Den(3)+tempa+tempb-tempc;
    Denex(7) = alpha*Den(7)+tempb+tempc;
    tempa = 2*Prod(knew,1)*Prod(knew,4);
    Denex(4) = alpha*Den(4)+tempa+Par(2)-Par(3);
    tempa = 2*Prod(knew,1)*Prod(knew,5);
    Denex(5) = alpha*Den(5)+tempa+Prod(knew,2)*Prod(knew,3);
    Denex(8) = alpha*Den(8)+Par(4)-Par(5);
    Denex(9) = alpha*Den(9)+Prod(knew,4)*Prod(knew,5);

    % Seek the value of the angle that maximizes the modulus of DENOM.
    tsum = Denex(1)+Denex(2)+Denex(4)+Denex(6)+Denex(8);
    denold = tsum;
    denmax = tsum;
    isave = 0;
    iu = 49;
    temp = 2*pi/(iu+1);
    Par(1) = 1;
    for i=1:iu
        angle = i*temp;
        Par(2) = cos(angle);
        Par(3) = sin(angle);
        for j=4:2:8
            Par(j) = Par(2)*Par(j-2)-Par(3)*Par(j-1);
            Par(j+1) = Par(2)*Par(j-1)+Par(3)*Par(j-2);
        end
        sumold = tsum;
        tsum = Denex(1:9)*Par(1:9)';
        if (abs(tsum)>abs(denmax))
            denmax = tsum;
            isave = i;
            tempa = sumold;
        elseif (i==isave+1)
            tempb = tsum;
        end
    end
    if (isave==0)
        tempa = tsum;
    end
    if (isave==iu)
        tempb = denold;
    end
    step = 0;
    if (tempa~=tempb)
        tempa = tempa-denmax;
        tempb = tempb-denmax;
        step = .5*(tempa-tempb)/(tempa+tempb);
    end
    angle = temp*(isave+step);

    % Calculate the new parameters of the denominator, the new Vlag vector
    % and the new D. Then test for convergence.
    Par(2) = cos(angle);
    Par(3) = sin(angle);
    for j=4:2:8
        Par(j) = Par(2)*Par(j-2)-Par(3)*Par(j-1);
        Par(j+1) = Par(2)*Par(j-1)+Par(3)*Par(j-2);
    end
    beta = Den(1:9)*Par(1:9)';
    denmax = Denex(1:9)*Par(1:9)';
    Vlag = Par(1:5)*Prod(1:npt+n,1:5)';
    tau = Vlag(knew);
    D(1:n) = Par(2)*D(1:n)+Par(3)*S(1:n);
    W(1:n) = Xopt(1:n)+D(1:n);
    dd = sum(D(1:n).^2);
    tempa = D(1:n)*W(1:n)';
    tempb = sum(W(1:n).^2);
    if (iterc>=n)
        break; % GOTO 340
    end
    if (iterc>1)
        densav = max(densav,denold);
    end
    if (abs(denmax)<=1.1*abs(densav))
        break; % GOTO 340
    end
    densav = denmax;

    % Set S to half the gradient of the denominator with respect to D.
    % Then branch for the next iteration.
    for i=1:n
        temp = tempa*Xopt(i)+tempb*D(i)-Vlag(npt+i);
        S(i) = tau*Bmat(knew,i)+alpha*temp;
    end
    %---------------------
    S = S+((tau*W(n+[1:npt])-alpha*Vlag(1:npt)).*(W(1:n)*Xpt'))*Xpt;
%     for k=1:npt
%         tsum = Xpt(k,1:n)*W(1:n)';
%         temp = (tau*W(n+k)-alpha*Vlag(k))*tsum;
%         S(1:n) = S(1:n)+temp*Xpt(k,1:n);
%     end
    %---------------------------
    ss = sum(S.^2);
    ds = D(1:n)*S(1:n)';
    ssden = dd*ss-ds*ds;
    if (ssden<(1e-8)*dd*ss)
        break; % GOTO 340
    end
end

% Set the vector W before the RETURN from the subroutine.
W(1:npt+n) = Wvec(1:npt+n,1:5)*Par(1:5)';
Vlag(kopt) = Vlag(kopt)+1;
return;
