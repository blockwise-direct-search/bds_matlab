function  [Step,crvmin] = trsapp(n,npt,Xopt,Xpt,Gq,Hq,Pq,delta)
% n is the number of variables of a quadratic objective function, Q say.
% The arguments npt, Xopt, Xpt, Gq, Hq and Pq have their usual meanings,
%   in order to define the current quadratic model Q.
% delta is the trust region radius, and has to be positive.
% Step will be set to the calculated trial step.
% crvmin will be set to the least curvature of H along the conjugate
%   directions that occur, except that it is set to zero if Step goes
%   all the way to the trust region boundary.
%
% The calculation of Step begins with the truncated conjugate gradient
% method. If the boundary of the trust region is reached, then further
% changes to Step may be made, each one being in the 2D space spanned
% by the current Step and the corresponding gradient of Q. Thus Step
% should provide a substantial reduction to Q within the trust region.

% Initialization, which includes setting Hd to H times Xopt.
delsq = delta*delta;
iterc = 0;
itermax = n;
itersw = itermax;
D(1:n) = Xopt(1:n);
Hd = hdset(n,npt,Pq,Xpt,D,Hq);

% Prepare for the first line search.
qred = 0;
Step = zeros(1,n);
Hs = zeros(1,n);
G = Gq+Hd;
D = -G;
dd = sum(D.^2);
crvmin = 0;
if (dd==0)
    return;
end
ds = 0;
ss = 0;
gg = dd;
ggbeg = gg;

% Calculate the step to the trust region boundary and the product Hd.
iterc = iterc+1;
temp = delsq-ss;
bstep = temp/(ds+sqrt(ds*ds+dd*temp));
Hd = hdset(n,npt,Pq,Xpt,D,Hq);

while (iterc<=itersw)
    flag40 = 1;
    dhd = D*Hd';

    % Update crvmin and set the step-length alpha.
    alpha = bstep;
    if (dhd>0)
        temp = dhd/dd;
        if (iterc==1)
            crvmin = temp;
        end
        crvmin = min(crvmin,temp);
        alpha = min(alpha,gg/dhd);
    end
    qadd = alpha*(gg-.5*alpha*dhd);
    qred = qred+qadd;

    % Update Step and Hs.
    ggsav = gg;
    Step = Step+alpha*D;
    Hs = Hs+alpha*Hd;
    gg = sum((G+Hs).^2);

    % Begin another conjugate direction iteration if required.
    if (alpha<bstep)
        if (qadd<=.01*qred)
            return;
        elseif (gg<=(1e-4)*ggbeg)
            return;
        elseif (iterc==itermax)
            return;
        end
        D = (gg/ggsav)*D-G-Hs;
        dd = sum(D.^2);
        ds = sum(D*Step');
        ss = sum(Step.^2);
        flag40 = 1;
        if (ds<=0)
            return;
        elseif (ss<delsq)
            flag40 = 0; % GOTO 40
            % Calculate the step to the trust region boundary and the product Hd.
            iterc = iterc+1;
            temp = delsq-ss;
            bstep = temp/(ds+sqrt(ds*ds+dd*temp));
            Hd = hdset(n,npt,Pq,Xpt,D,Hq);
        end
    end

    if (flag40)
        crvmin = 0;
        itersw = iterc;

        % Test whether an alternative iteration is required.
        if (gg<=(1e-4)*ggbeg)
            return;
        else
            sg = Step*G';
            shs = Step*Hs';
            sgk = sg+shs;
            angtest = sgk/sqrt(gg*delsq);
            if (angtest<=-.99)
                return;
            end
        end
        % Begin the alternative iteration by calculating D and Hd and some
        % scalar products.
        iterc = iterc+1;
        temp = sqrt(delsq*gg-sgk*sgk);
        D = (delsq/temp)*(G+Hs)-(sgk/temp)*Step;
        Hd = hdset(n,npt,Pq,Xpt,D,Hq);
    end
end % end of first while loop

while (1)
    dg = D*G';
    dhd = Hd*D';
    dhs = Hd*Step';

    % Seek the value of the angle that minimizes Q.
    cf = .5*(shs-dhd);
    qbeg = sg+cf;
    qsav = qbeg;
    qmin = qbeg;
    isave = 0;
    iu = 49;
    temp = 2*pi/(iu+1);
    for i=1:iu
        angle = i*temp;
        cth = cos(angle);
        sth = sin(angle);
        qnew = (sg+cf*cth)*cth+(dg+dhs*cth)*sth;
        if (qnew<qmin)
            qmin = qnew;
            isave = i;
            tempa = qsav;
        elseif (i==isave+1)
            tempb = qnew;
        end
        qsav = qnew;
    end
    if (isave==0)
        tempa = qnew;
    end
    if (isave==iu)
        tempb = qbeg;
    end
    angle = 0;
    if (tempa~=tempb)
        tempa = tempa-qmin;
        tempb = tempb-qmin;
        angle = .5*(tempa-tempb)/(tempa+tempb);
    end
    angle = temp*(isave+angle);

    % Calculate the new Step and Hs. Then test for convergence.
    cth = cos(angle);
    sth = sin(angle);
    reduc = qbeg-(sg+cf*cth)*cth-(dg+dhs*cth)*sth;
    Step = cth*Step+sth*D;
    Hs = cth*Hs+sth*Hd;
    gg = sum((G+Hs).^2);
    qred = qred+reduc;
    ratio = reduc/qred;
    if (iterc<itermax) && (ratio>.01)
        % Test whether an alternative iteration is required.
        if (gg<=(1e-4)*ggbeg)
            return;
        else
            sg = Step*G';
            shs = Step*Hs';
            sgk = sg+shs;
            angtest = sgk/sqrt(gg*delsq);
            if (angtest<=-.99)
                return;
            end
        end
        % Begin the alternative iteration by calculating D and Hd and some
        % scalar products.
        iterc = iterc+1;
        temp = sqrt(delsq*gg-sgk*sgk);
        D = (delsq/temp)*(G+Hs)-(sgk/temp)*Step;
        Hd = hdset(n,npt,Pq,Xpt,D,Hq);
    else
        return;
    end
end % end of second while loop

    function Hd = hdset(n,npt,Pq,Xpt,D,Hq)
        % The following instructions act as a subroutine for setting the vector
        % Hd to the vector D multiplied by the second derivative matrix of Q.
        % They are called from three different places, which are distinguished
        % by the value of iterc.
        
        Hd = ((D*Xpt').*Pq)*Xpt;
        %         Hd = zeros(1,n);
        %         for k=1:npt
        %             temp = Pq(k)*(Xpt(k,1:n)*D(1:n)');
        %             Hd(1:n) = Hd(1:n)+temp*Xpt(k,1:n);
        %         end

        Hd(1) = Hd(1) + D(1:n)*Hq(1+.5*(1:n).*((1:n)-1))';
        for k=2:n
           Hd(k) = Hd(k)+D(1:k-1)*Hq((1:k-1)+.5*k*(k-1))'+D(k:n)*Hq(k+.5*(k:n).*((k:n)-1))';
        end
%         
%         ih = 0;
%         for j=1:n
%             for i=1:j
%                 ih = ih+1;
%                 if (i<j)
%                      Hd(j) = Hd(j)+Hq(ih)*D(i);
%                 end
%                 Hd(i) = Hd(i)+Hq(ih)*D(j);
%             end
%         end
        %norm((Hd-Hd1)./abs(Hd))
    end % End of hdset subroutine

end % end of trsapp function
