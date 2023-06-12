% sigma = eye(5);
% for i = 1:5
%     sigma(i,i) = 20^i;
% end
%sigma(5,5) = 10^7;
%p = macup('PALMER2C');
sample = [-1.745329000000000
  -1.570796000000000
  -1.396263000000000
  -1.221730000000000
  -1.047198000000000
  -0.937187000000000
  -0.872665000000000
  -0.698132000000000
  -0.523599000000000
  -0.349066000000000
  -0.174533000000000
                   0
   0.174533000000000
   0.349066000000000
   0.523599000000000
   0.698132000000000
   0.872665000000000
   0.937187000000000
   1.047198000000000
   1.221730000000000
   1.396263000000000
   1.570796000000000
   1.745329000000000];
%X = [ones(23,1) sample.^2 sample.^4 sample.^6 sample.^8 sample.^10 sample.^12 sample.^14];
X = [ones(23,1) sample.^2 sample.^4];
XX = X'*X;
e = eig(XX);
n = length(e);
cond = e(n)/e(1)
[VEC,D] = eig(XX)
%e(8) = 10*e(8);
V = diag(e);
%rng(100)
QA = randn(n,n);
[Q, R] = qr(QA);
%rng(100)
x0 = randn(n,1);
%norm(x0)
fun = @(x) x'*XX*x;
options.nb = 1;
options.maxfun = 2e4;
options.polling_inner = "complete";
%fun(x0)
options.shrink = 0.5;
[xval, fval, exitflag, output] = blockwise_direct_search(fun, x0, options)
com_fhist = output.fhist;
%xval
options.expand = 2;
options.shrink = 0.5;
options.polling_inner = "opportunistic";
[xval, fval, exitflag, output] = blockwise_direct_search(fun, x0, options)
%xval
opp_fhist = output.fhist;
min_length = min([length(com_fhist), 5000]);
semilogy(com_fhist(1:min_length), 'b');
hold on
min_length = min([length(opp_fhist), 5000]);
semilogy(opp_fhist(1:min_length), 'r')
legend('com', 'opp')
