%% Test files 
%%1-D poisson with sin(k*pi*x) right hand side
n=100;
x = linspace(0,1,n+2);
x = x(2:end-1);

A = zeros(n,n);
A(1:(n+1):n^2) = 2;
A(2:(n+1):n^2) = -1;
A((n+1):(n+1):n^2) = -1;

A = A + eye(n);

b = sin(10.3*pi*x);
b = b(:);

Meurant(A,b,1:(ceil(n/10)):n);


%% 2- D poisson with all ones right hand side
n=144;
x = linspace(0,1,n+2);
x = x(2:end-1);

A = gallery('poisson',sqrt(n));

A = A + 0.5*eye(n);
A = A;

b = ones(n,1);
b = b(:);

Meurant(A,b,1:1:30);

%% random squared example
n=20;
A = rand(n);
A = A'*A;

A = eye(n) + A;
b = rand(n,1);

Meurant(A,b,1:10:n);


%% Chebyshev diff example
%% You need chebfun to run example

% 
% %%Define projection
% P = legpoly(0:n-1);
% x = chebfun('x'); 
% [Q, R] = qr( repmat(1-x.^2,1,n).*P  ); 
% 
% Qc = chebcoeffs(Q);
% P = Qc*((Qc(1:end-2,:))');
% 
% D = ChebDiffMat(n+2);
% L = - D*D;
% %L = L(1:end-2,:);
% 
% b = chebfun('exp(sin(1*pi*x))',size(chebfun,);
% 
% bc = chebcoeffs(b,size(Qc,1));
% 
% norm(chebcoeffs(Q*(Q'*b)) - Qc*(Qc'*bc))
% 

