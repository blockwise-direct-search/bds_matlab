% Mnewuoa.m Created 6/5/08 by Stefan Wild and Jorge More' (updated 1/5/12)
% 
%
% This is Matlab code transcribed from the original NEWUOA fortran code of 
% M.J.D. Powell (mjdp@cam.ac.uk) [dated December 16th, 2004], please see 
% the included README. No warranty is provided, please do not distribute
% this code without emailing wild@mcs.anl.gov.
%
% INPUTS:  ----------------------------------------------------------------
%  n must be set to the number of variables and must be at least two.
%  npt is the number of interpolation conditions. Its value must be in the
%    interval [n+2,(n+1)(n+2)/2].
%  Initial values of the variables must be set in X(1),X(2),...,X(n). They
%   will be changed to the values that give the least calculated f.
% rhobeg and rhoend must be set to the initial and final values of a trust
%   region radius, so both must be positive with rhoend<=rhobeg. Typically
%   rhobeg should be about one tenth of the greatest expected change to a
%   variable, and rhoend should indicate the accuracy that is required in
%   the final values of the variables.
%  The value of iprint should be set to 0, 1, 2 or 3, which controls the
%    amount of printing. Specifically, there is no output if iprint=0 and
%   there is output only at the return if iprint=1. Otherwise, each new
%   value of rho is printed, with the best vector of variables so far and
%    the corresponding value of the objective function. Further, each new
%    value of f with its variables are output if iprint=3.
%  maxfun must be set to an upper bound on the number of calls of calfun.
%  SUBROUTINE f=calfun(X) must be provided by the user. It must set f to
% the value of the objective function for the variables X(1),X(2),...,X(n).
function X=mnewuoa(calfun,n,npt,X,rhobeg,rhoend,iprint,maxfun)

if (npt<n+2) || (npt>(n+2)*(n+1)/2)
    disp('Return from NEWUOA because npt is not in the required interval')
    return
end

Xbase = X; % Assumed to be a row vector?
% Xbase will hold a shift of origin that should reduce the contributions
%   from rounding errors to values of the model and Lagrange functions.
Xopt = zeros(1,n);
% Xopt will be set to the displacement from Xbase of the vector of
%   variables that provides the least calculated f so far.
%Xnew = zeros(1,n);
% Xnew will be set to the displacement from Xbase of the vector of
%   variables for the current calculation of f.
Xpt = zeros(npt,n);
% Xpt will contain the interpolation point coordinates relative to Xbase.
Fval = zeros(1,npt);
% Fval will hold the values of f at the interpolation points.
Gq = zeros(1,n);
% Gq will hold the gradient of the quadratic model at Xbase.
Hq = zeros(1,n*(n+1)/2);
% Hq will hold the explicit second derivatives of the quadratic model.
Pq = zeros(1,npt); % SHould be a column vector
% Pq will contain the parameters of the implicit second derivatives of
%   the quadratic model.
Bmat = zeros(npt+n,n);
% Bmat will hold the last n columns of H.
Zmat = zeros(npt,npt-n-1);
% Zmat will hold the factorization of the leading npt by npt submatrix of
%   H, this factorization being Zmat times Diag(DZ) times Zmat^T, where
%   the elements of DZ are plus or minus one, as specified by idz.
%D = zeros(1,n);
% D is reserved for trial steps from Xopt.
Vlag = zeros(1,npt+n);
% Vlag will contain the values of the Lagrange functions at a new point X.
%   They are part of a product that requires Vlag to be of length ndim.
W = zeros(1,10*(npt+n));
% The array W will be used for working space.

% Begin the initialization procedure. nf becomes one more than the number
% of function values so far. The coordinates of the displacement of the
% next initial interpolation point from Xbase are set in Xpt(nf,.).
rhosq = rhobeg*rhobeg;
recip = 1/rhosq;
reciq = sqrt(.5)/rhosq;
nf = 0;
while nf<npt
    nfm = nf;
    nfmm = nf-n;
    nf = nf+1;
    if (nfm<=2*n)
        if (nfm>=1) && (nfm<=n)
            Xpt(nf,nfm) = rhobeg;
        elseif (nfm>n)
            Xpt(nf,nfmm) = -rhobeg;
        end
    else
        itemp = floor((nfmm-1)/n);
        jpt = nfm-itemp*n-n;
        ipt = jpt+itemp;
        if (ipt>n)
            itemp = jpt;
            jpt = ipt-n;
            ipt = itemp;
        end
        xipt = rhobeg;
        if (Fval(ipt+n+1)<Fval(ipt+1))
            xipt = -xipt;
        end
        xjpt = rhobeg;
        if (Fval(jpt+n+1)<Fval(jpt+1))
            xjpt = -xjpt;
        end
        Xpt(nf,ipt) = xipt;
        Xpt(nf,jpt) = xjpt;
    end
    % Calculate the next value of f. The least function value so far
    % and its index are required.
    X(1:n) = Xpt(nf,1:n) + Xbase(1:n);

    if (nf>max(1,maxfun))
        nf = nf-1;
        if (iprint>0)
            disp('Return from NEWUOA because');
            disp('calfun has been called maxfun times.');
        end
        if (fopt<= f)
            X(1:n) = Xbase(1:n)+Xopt(1:n);
            f = fopt;
        end
        if (iprint>=1)
            disp('At the return from NEWUOA');
            disp(strcat('    Number of function values =',num2str(nf)));
            disp(f)
            disp(X)
        end
        return;
    end
    f = calfun(X);
    if (iprint==3)
        disp(strcat('Function number',num2str(nf),'    f =',num2str(f)));
        disp('    The corresponding X is:');
        disp(X)
    end
    Fval(nf) = f;
    if (nf==1)
        fbeg = f;
        fopt = f;
        kopt = 1;
    elseif (f<fopt)
        fopt = f;
        kopt = nf;
    end

    % Set the nonzero initial elements of Bmat and the quadratic model in
    % the cases when nf is at most 2*n+1.
    if (nfm<=2*n)
        if (nfm>=1) && (nfm<=n)
            Gq(nfm) = (f-fbeg)/rhobeg;
            if (npt<nf+n)
                Bmat([1, nf],nfm) = [-1/rhobeg; 1/rhobeg];
                Bmat(npt+nfm,nfm) = -.5*rhosq;
            end
        elseif (nfm>n)
            Bmat([nf-n nf],nfmm) = .5/rhobeg*[1; -1];
            Zmat([1 nf-n nf],nfmm) = [-2*reciq; reciq; reciq];
            temp = (fbeg-f)/rhobeg;
            Hq((nfmm*(nfmm+1))/2) = (Gq(nfmm)-temp)/rhobeg;
            Gq(nfmm) = .5*(Gq(nfmm)+temp);
        end

        % Set the off-diagonal second derivatives of the Lagrange functions
        % and the initial quadratic model.
    else
        ih = (ipt*(ipt-1))/2+jpt;
        if (xipt<0)
            ipt = ipt+n;
        end
        if (xjpt<0)
            jpt = jpt+n;
        end
        Zmat([1, nf, ipt+1, jpt+1],nfmm) = recip*[1; 1; -1; -1];
        Hq(ih) = (fbeg-Fval(ipt+1)-Fval(jpt+1)+f)/(xipt*xjpt);
    end
end

% Begin the iterative procedure, because the initial model is complete.
rho = rhobeg;
delta = rho;
idz = 1;
diffa = 0;
diffb = 0;
itest = 0;
Xopt(1:n) = Xpt(kopt,1:n);
xoptsq = sum(Xopt.^2);
nfsav = nf;

loop100 = 1; % Indicator to proceed to loop100
while (1)
    % Generate the next trust region step and test its length. Set knew
    % to -1 if the purpose of the next f will be to improve the model.
    if (loop100)
        loop120 = 1; % Indicator to proceed to loop120
        knew = 0;
        [D,crvmin] = trsapp(n,npt,Xopt,Xpt,Gq,Hq,Pq,delta);
        dsq = sum(D(1:n).^2);
        dnorm = min(delta,sqrt(dsq));
        if (dnorm<.5*rho)
            knew = -1;
            delta = .1*delta;
            ratio = -1;
            if (delta<=1.5*rho)
                delta = rho;
            end

            temp = 0.125*crvmin*rho*rho;
            if (nf<=nfsav+2)
                loop120 = 0; loop290 = 0; loop460 = 1; % GOTO 460
            elseif (temp <= max([diffa,diffb,diffc]))
                loop120 = 0; loop290 = 0; loop460 = 1; % GOTO 460
            elseif (rho>rhoend)
                delta = .5*rho;
                ratio = rho/rhoend;
                if (ratio<=16)
                    rho = rhoend;
                elseif (ratio<=250)
                    rho = sqrt(ratio)*rhoend;
                else
                    rho = .1*rho;
                end
                delta = max(delta,rho);
                if (iprint>=2)
                    disp(strcat('New rho =',num2str(rho)));
                    disp(strcat('Number of function values=',num2str(nf)));
                    disp(strcat('Least value of f =',num2str(f)));
                    disp('The corresponding X is:');
                    disp(Xbase+Xopt);
                end
                nfsav = nf;
                loop120 = 0; loop290 = 0; loop460 =0; % Return to loop100
            end
        end
    end % End of loop 100

    if (loop120)
        % Shift Xbase if Xopt may be too far from Xbase. First make the
        % changes to Bmat that do not depend on Zmat.
        if (dsq<=.001*xoptsq)
            tempq = .25*xoptsq;
            for k=1:npt
                tsum = Xpt(k,1:n)*Xopt(1:n)';
                temp = Pq(k)*tsum;
                tsum = tsum-.5*xoptsq;
                W(npt+k) = tsum;
                for i=1:n
                    Gq(i) = Gq(i)+temp*Xpt(k,i);
                    Xpt(k,i) = Xpt(k,i)-.5*Xopt(i);
                    Vlag(i) = Bmat(k,i);
                    W(i) = tsum*Xpt(k,i)+tempq*Xopt(i);
                    ip = npt+i;
                    Bmat(ip,1:i) = Bmat(ip,1:i)+Vlag(i)*W(1:i)+W(i)*Vlag(1:i);
                end
            end

            % Then the revisions of Bmat that depend on Zmat are calculated.
            for k=1:(npt-n-1)
                W(1:npt) = W(npt+[1:npt])*Zmat(1:npt,k);
                sumz = sum(Zmat(1:npt,k));
                for j=1:n
                    tsum = tempq*sumz*Xopt(j);
                    tsum = tsum + W(1:npt)*Xpt(1:npt,j);
                    Vlag(j) = tsum;
                    if (k<idz)
                        tsum = -tsum;
                    end
                    Bmat(1:npt,j) = Bmat(1:npt,j)+tsum*Zmat(1:npt,k);
                end
                for i=1:n
                    ip = i+npt;
                    temp = Vlag(i);
                    if  (k<idz)
                        temp = -temp;
                    end
                    Bmat(ip,1:i) = Bmat(ip,1:i)+temp*Vlag(1:i);
                end
            end

            % The following instructions complete the shift of Xbase,
            % including the changes to the parameters of the model.
            ih = 0;
            for j=1:n
                W(j) = Pq(1:npt)*Xpt(1:npt,j);
                Xpt(1:npt,j) = Xpt(1:npt,j)-.5*Xopt(j);
                for i=1:j
                    ih = ih+1;
                    if (i<j)
                        Gq(j) = Gq(j)+Hq(ih)*Xopt(i);
                    end
                    Gq(i) = Gq(i)+Hq(ih)*Xopt(j);
                    Hq(ih) = Hq(ih)+W(i)*Xopt(j)+Xopt(i)*W(j);
                    Bmat(npt+i,j) = Bmat(npt+j,i);
                end
            end
            Xbase(1:n) = Xbase(1:n)+Xopt(1:n);
            Xopt = zeros(1,n);
            xoptsq = 0;
        end

        % Pick the model step if knew is positive. A different choice of D
        % may be made later, if the choice of D by BIGLAG causes
        % substantial cancellation in DENOM.
        if (knew>0)
            [D,alpha] = biglag(n,npt,Xopt,Xpt,Bmat,Zmat,idz,knew,delta);
        end

        % Calculate Vlag and beta for the current choice of D. The first
        % npt components of W_check will be held in W.
        
        
%         %-----------------------
       W(1:npt) = (D*Xpt').*((.5*D+Xopt)*Xpt');
       Vlag(1:npt) = D*Bmat(1:npt,1:n)';
%         for k=1:npt
%             suma = Xpt(k,1:n)*D(1:n)';
%             sumb = Xpt(k,1:n)*Xopt(1:n)';
%             tsum = Bmat(k,1:n)*D(1:n)';
%             W(k) = suma*(.5*suma+sumb);
%             Vlag(k) = tsum;
%         end
        %-----------------------
        
        
        
        %-----------------------
        beta = W(1:npt)*Zmat(1:npt,1:idz-1)*Zmat(1:npt,1:idz-1)'*W(1:npt)';
        Vlag(1:npt) = Vlag(1:npt)-W(1:npt)*Zmat(1:npt,1:idz-1)*Zmat(1:npt,1:idz-1)';
        beta = beta-W(1:npt)*Zmat(1:npt,idz:npt-n-1)*Zmat(1:npt,idz:npt-n-1)'*W(1:npt)';
        tsum = W(1:npt)*Zmat(1:npt,idz:npt-n-1);
        Vlag(1:npt) = Vlag(1:npt)+tsum*Zmat(1:npt,idz:npt-n-1)';

%         beta = 0;
%         for k=1:(npt-n-1)
%             tsum = W(1:npt)*Zmat(1:npt,k);
%             if (k<idz)
%                 beta = beta+tsum*tsum;
%                 tsum = -tsum;
%             else
%                 beta = beta-tsum*tsum;
%             end
%             Vlag(1:npt) = Vlag(1:npt)+tsum*Zmat(1:npt,k)';
%         end
        %-----------------------
        
        %-----------------
%         bsum = 0;
%         dx = 0;
%         for j=1:n
%             tsum = W(1:npt)*Bmat(1:npt,j);
%             bsum = bsum+tsum*D(j);
%             jp = npt+j;
%             tsum = tsum+Bmat(jp,1:n)*D(1:n)';
%             Vlag(jp) = tsum;
%             bsum = bsum+tsum*D(j);
%             dx = dx+D(j)*Xopt(j);
%         end
            bsum = W(1:npt)*Bmat(1:npt,1:n)*D(1:n)';
            tsum = W(1:npt)*Bmat(1:npt,1:n)+D(1:n)*Bmat(npt+(1:n),1:n)';
            Vlag(npt+(1:n)) = tsum;
            bsum = bsum+tsum*D(1:n)';
            dx = D(1:n)*Xopt(1:n)';       
            %-----------------
        
        beta = dx*dx+dsq*(xoptsq+dx+dx+.5*dsq)+beta-bsum;
        Vlag(kopt) = Vlag(kopt)+1;

        % If knew is positive and if the cancellation in DENOM is
        % unacceptable, then BIGDEN calculates an alternative model step,
        % Xnew being used for working space.
        if (knew>0)
            temp = 1+alpha*beta/Vlag(knew)^2;
            if (abs(temp)<=0.8)
                [D,W,Vlag,beta] = bigden(n,npt,Xopt,Xpt,Bmat,Zmat,idz,kopt,knew,D);
            end
        end
        loop290 = 1; % Indicator to proceed to loop290
    end % End of loop 120

    if (loop290)
        %  Calculate the next value of the objective function.
        Xnew(1:n) = Xopt(1:n)+D(1:n);
        X(1:n) = Xbase(1:n)+Xnew(1:n);
        nf = nf+1;

        if (nf>max(1,maxfun))
            nf = nf-1;
            if (iprint>0)
                disp('Return from NEWUOA because ');
                disp('calfun has been called maxfun times.');
            end
            if (fopt<= f)
                X(1:n) = Xbase(1:n)+Xopt(1:n);
                f = fopt;
            end
            if (iprint>=1)
                disp('At the return from NEWUOA');
                disp(strcat('   Number of function values =',num2str(nf)));
                disp(f)
                disp(X)
            end
            return;
        end
        f = calfun(X);
        if (iprint==3)
            disp(strcat('Function number',num2str(nf),'  f =',num2str(f)));
            disp('    The corresponding X is:');
            disp(X)
        end

        if (knew==-1)
            if (fopt<=f)
                X(1:n) = Xbase(1:n)+Xopt(1:n);
                f = fopt;
            end
            if (iprint>=1)
                disp(strcat('At the return from NEWUOA'));
                disp(strcat('Number of function values =',num2str(nf)));
                disp(f)
                disp(X)
            end
            return;
        end

        % Use the quadratic model to predict the change in f due to the
        % step D, and set diff to the error of this prediction.
        vquad = 0;
        ih = 0;
        for j=1:n
            vquad = vquad+D(j)*Gq(j);
            for i=1:j
                ih = ih+1;
                temp = D(i)*Xnew(j)+D(j)*Xopt(i);
                if (i==j)
                    temp = .5*temp;
                end
                vquad = vquad+temp*Hq(ih);
            end
        end
        vquad = vquad+sum(Pq(1:npt).*W(1:npt));

        diff = f-fopt-vquad;
        diffc = diffb;
        diffb = diffa;
        diffa = abs(diff);
        if (dnorm>rho)
            nfsav = nf;
        end

        % Update fopt and Xopt if the new f is the least value of the
        % objective function so far. The branch when knew is positive
        % occurs if D is not a trust region step.
        fsave = fopt;
        if (f<fopt)
            fopt=f;
            Xopt(1:n) = Xnew(1:n);
            xoptsq = sum(Xopt.^2);
        end
        ksave = knew;
        if (knew<=0)
            % Pick the next value of delta after a trust region step.
            if (vquad>=0)
                if (iprint>0)
                    disp('Return from NEWUOA because a trust');
                    disp(' region step has failed to reduce Q.');
                end
                if (fopt<=f)
                    X(1:n) = Xbase(1:n)+Xopt(1:n);
                    f = fopt;
                end
                if (iprint>=1)
                    disp(strcat('At the return from NEWUOA'));
                    disp(strcat('Number of function values=',num2str(nf)));
                    disp(f)
                    disp(X)
                end
                return;
            end
            ratio = (f-fsave)/vquad;
            if (ratio<=.1)
                delta = .5*dnorm;
            elseif (ratio<=.7)
                delta = max(.5*delta,dnorm);
            else
                delta = max(.5*delta,dnorm+dnorm);
            end
            if (delta<=1.5*rho)
                delta = rho;
            end

            % Set knew to the index of the next interpolation point to be
            % deleted.
            rhosq = max(.1*delta,rho)^2;
            ktemp = 0;
            detrat = 0;
            if (f>=fsave)
                ktemp = kopt;
                detrat = 1;
            end
            for k=1:npt
                hdiag = 0;
                for j=1:(npt-n-1)
                    temp = 1;
                    if (j<idz)
                        temp = -1;
                    end
                    hdiag = hdiag+temp*Zmat(k,j)^2;
                end
                temp = abs(beta*hdiag+Vlag(k)^2);
                distsq = sum((Xpt(k,1:n)-Xopt(1:n)).^2);
                if (distsq>rhosq)
                    temp = temp*(distsq/rhosq)^3;
                end
                if (temp>detrat) && (k~=ktemp)
                    detrat = temp;
                    knew = k;
                end
            end
        end
        if (knew==0)
            loop460 = 1; % GOTO 460
        else
            % Update Bmat, Zmat and idz, so that the knew-th interpolation
            % point can be moved. Begin the updating of the quadratic
            % model, starting with the explicit second derivative term.
            [Bmat,Zmat,idz,knew] = update_models(n,npt,Bmat,Zmat,idz,Vlag,beta,knew);
            Fval(knew) = f;
            ih = 0;
            for i=1:n
                temp = Pq(knew)*Xpt(knew,i);
                for j=1:i
                    ih = ih+1;
                    Hq(ih) = Hq(ih)+temp*Xpt(knew,j);
                end
            end
            Pq(knew) = 0;

            % Update the other 2nd derivative parameters, then the gradient
            % vector of the model. Also include the new interpolation point.
            %--------------------------
%             for j=1:(npt-n-1)
%                 temp = diff*Zmat(knew,j);
%                 if (j<idz)
%                     temp = -temp;
%                 end
%                 Pq(1:npt) = Pq(1:npt)+temp*Zmat(1:npt,j)';
%             end
            Pq(1:npt) = Pq(1:npt)-diff*Zmat(knew,1:idz-1)*Zmat(1:npt,1:idz-1)';
            Pq(1:npt) = Pq(1:npt)+diff*Zmat(knew,idz:npt-n-1)*Zmat(1:npt,idz:npt-n-1)';
            %--------------------------
            Gq(1:n) = Gq(1:n)+diff*Bmat(knew,1:n);
            gqsq = sum(Gq(1:n).^2);
            Xpt(knew,1:n) = Xnew(1:n);

            % If a trust region step makes a small change to the objective
            % function, then calculate the gradient of the least Frobenius
            % norm interpolant at Xbase, and store it in W, using Vlag for
            % a vector of right hand sides.
            if (ksave==0) && (delta==rho)
                if (abs(ratio)>.01)
                    itest = 0;
                else
                    Vlag(1:npt) = Fval(1:npt)-Fval(kopt);
                    gisq = 0;
                    for i=1:n
                        tsum = Vlag(1:npt)*Bmat(1:npt,i);
                        gisq = gisq+tsum*tsum;
                        W(i) = tsum;
                    end

                    % Test whether to replace the new quadratic model by
                    % the least Frobenius norm interpolant, making the
                    % replacement if the test is satisfied.
                    itest = itest+1;
                    if (gqsq<100*gisq)
                        itest = 0;
                    end
                    if (itest>=3)
                        Gq(1:n) = W(1:n);
                        Hq = zeros(1,(n*(n+1))/2);
                        for j=1:(npt-n-1)
                            W(j) = Vlag(1:npt)*Zmat(1:npt,j);
                            if (j<idz)
                                W(j) = -W(j);
                            end
                        end
                        for k=1:npt
                            Pq(k) = Zmat(k,1:(npt-n-1))*W(1:(npt-n-1))';
                        end
                        itest = 0;
                    end
                end
            end
            if (f<fsave)
                kopt = knew;
            end

            % If a trust region step has provided a sufficient decrease in
            % f, then branch for another trust region calculation. The case
            % ksave>0 occurs when the new function value was calculated by
            % a model step.
            if (f <= fsave+.1*vquad)
                loop460 = 0; loop100 = 1; % GOTO 100
            elseif (ksave>0)
                loop460 = 0; loop100 = 1; % GOTO 100
            else
                loop460 = 1; % continue on to loop 460
            end

            % Alternatively, find out if the interpolation points are close
            % enough to the best point so far.
            knew = 0;
        end
    end % End of loop290

    if (loop460)
        distsq=4*delta*delta;
        for k=1:npt
            tsum = sum((Xpt(k,1:n)-Xopt(1:n)).^2);
            if (tsum>distsq)
                knew = k;
                distsq = tsum;
            end
        end

        % If knew is positive, then set dstep, and branch back for the next
        % iteration, which will generate a "model step".
        if (knew>0)
            dstep = max(min(.1*sqrt(distsq),.5*delta),rho);
            dsq = dstep*dstep;
            loop100 = 0; loop120 = 1; % GOTO 120
        elseif (ratio>0)
            loop100 = 1; % GOTO 100
        elseif (max(delta,dnorm)>rho)
            loop100 = 1; % GOTO 100

            % The calculations with the current value of rho are complete.
            % Pick the next values of rho and delta.
        elseif (rho>rhoend)
            delta = .5*rho;
            ratio = rho/rhoend;
            if (ratio<=16)
                rho = rhoend;
            elseif (ratio<=250)
                rho = sqrt(ratio)*rhoend;
            else
                rho = .1*rho;
            end
            delta = max(delta,rho);
            if (iprint>=2)
                disp(strcat('New rho =',num2str(rho)));
                disp(strcat('Number of function values =',num2str(nf)));
                disp(strcat('Least value of f =',num2str(f)))
                disp('The corresponding X is:')
                disp(Xbase+Xopt)
            end
            nfsav = nf;
            loop100 = 1; % GOTO 100

            % Return from the calculation, after another Newton-Raphson
            % step, if it is too short to have been tried before.
        elseif (knew==-1)
            loop100 = 0; loop120 = 0; loop290 = 1; % GOTO 290
        else
            if (fopt<=f)
                X(1:n) = Xbase(1:n)+Xopt(1:n);
                f = fopt;
            end
            if (iprint>=1)
                disp(strcat('At the return from NEWUOA'));
                disp(strcat(' Number of function values = ',num2str(nf)));
                disp(strcat(' f = ',num2str(f)));
                disp(strcat(' X = ',num2str(X)));
            end
            return;
        end
    end % End of loop460
end % End of while loop
