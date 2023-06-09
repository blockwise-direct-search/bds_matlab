function [Bmat,Zmat,idz,knew] = update_models(n,npt,Bmat,Zmat,idz,Vlag,beta,knew)
% The arrays Bmat and Zmat with idz are updated, in order to shift the
% interpolation point that has index knew. On entry, Vlag contains the
% components of the vector Theta*Wcheck+e_b of the updating formula
% (6.11), and beta holds the value of the parameter that has this name.

% Apply the rotations that put zeros in the knew-th row of Zmat.
jl = 1;
for j=2:npt-n-1
    if (j==idz)
        jl = idz;
    elseif (Zmat(knew,j)~=0)
        temp = sqrt(Zmat(knew,jl)^2+Zmat(knew,j)^2);
        tempa = Zmat(knew,jl)/temp;
        tempb = Zmat(knew,j)/temp;

%         for i=1:npt
%             temp = tempa*Zmat(i,jl)+tempb*Zmat(i,j);
%             Zmat(i,j) = tempa*Zmat(i,j)-tempb*Zmat(i,jl);
%             Zmat(i,jl) = temp;
%         end
        T = tempa*Zmat(:,jl)+tempb*Zmat(:,j);
        Zmat(:,j) = tempa*Zmat(:,j)-tempb*Zmat(:,jl);
        Zmat(:,jl) = T;
        
        Zmat(knew,j) = 0;
    end
end

% Put the first npt components of the knew-th column of HLAG into W,
% and calculate the parameters of the updating formula.
tempa = Zmat(knew,1);
if (idz>=2)
    tempa = -tempa;
end
if (jl>1)
    tempb = Zmat(knew,jl);
end
W = tempa*Zmat(:,1)';
if (jl>1)
    W = W+tempb*Zmat(:,jl)';
end
alpha = W(knew);
tau = Vlag(knew);
tausq = tau*tau;
denom = alpha*beta+tausq;
Vlag(knew) = Vlag(knew)-1;

% Complete the updating of Zmat when there is only one nonzero element
% in the knew-th row of the new matrix Zmat, but, if iflag is set to one,
% then the first column of Zmat will be exchanged with another one later.
iflag = 0;
if (jl==1)
    temp = sqrt(abs(denom));
    tempb = tempa/temp;
    tempa = tau/temp;
    Zmat(1:npt,1) = tempa*Zmat(1:npt,1)-tempb*Vlag(1:npt)';
    if (idz==1) && (temp<0)
        idz = 2;
    end
    if (idz>=2) && (temp>=0)
        iflag = 1;
    end
else

    % Complete the updating of Zmat in the alternative case.
    ja = 1;
    if (beta>=0)
        ja = jl;
    end
    jb = jl+1-ja;
    temp = Zmat(knew,jb)/denom;
    tempa = temp*beta;
    tempb = temp*tau;
    temp = Zmat(knew,ja);
    scala = 1/sqrt(abs(beta)*temp*temp+tausq);
    scalb = scala*sqrt(abs(denom));
    Zmat(1:npt,ja) = scala*(tau*Zmat(1:npt,ja)-temp*Vlag');
    Zmat(1:npt,jb) = scalb*(Zmat(1:npt,jb)-tempa*W'-tempb*Vlag');
    if (denom<=0)
        if (beta<0)
            idz = idz+1;
        else
            iflag = 1;
        end
    end
end

% idz is reduced in the following case, and usually the first column
% of Zmat is exchanged with a later one.
if (iflag==1)
    idz = idz-1;
    for i=1:npt
        temp = Zmat(i,1);
        Zmat(i,1) = Zmat(i,idz);
        Zmat(i,idz) = temp;
    end
end

% Finally, update the matrix Bmat.
for j=1:n
    jp = npt+j;
    W(jp) = Bmat(knew,j);
    tempa = (alpha*Vlag(jp)-tau*W(jp))/denom;
    tempb = (-beta*W(jp)-tau*Vlag(jp))/denom;

    %--------------
    for i=1:jp
        Bmat(i,j) = Bmat(i,j)+tempa*Vlag(i)+tempb*W(i);
        if (i>npt)
            Bmat(jp,i-npt) = Bmat(i,j);
        end
    end
%     Bmat(1:npt,j) = Bmat(1:npt,j)+tempa*Vlag(1:npt)'+tempb*W(1:npt)';
%     Bmat(npt+1:jp,j) = Bmat(npt+1:jp,j)+tempa*Vlag(npt+1:jp)'+tempb*W(npt+1:jp)';
%     Bmat(npt+j,1:j) = Bmat(npt+1:npt+j,j)';
    %--------------
end
return;
