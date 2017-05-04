%This is a script which compute the iterates of CG and the error in the 
% A norm at each step. Not meant to be efficient in any way
%
% cgerr2 : squared error at the kth step;
function [x ,cgerr] = VanillaCG(A,b)

    n = size(A,1);
    x0 = zeros(n,1);
    x(:,1) =x0;
    r(:,1) = b - A*x0;
    p(:,1) = r(:,1);

    sol = A\b;
    niter=n;
    
    cgerr(1) = sol'*A*sol;

    for k=1:niter
        alpha(k) = (r(:,k)'*r(:,k))/(p(:,k)'*A*p(:,k));
        x(:,k+1) = x(:,k)+alpha(k)*p(:,k);
        r(:,k+1) = r(:,k) - alpha(k)*A*p(:,k);
        beta(k) = (r(:,k+1)'*r(:,k+1))/(r(:,k)'*r(:,k));
        p(:,k+1) = r(:,k+1) + beta(k)*p(:,k);

        cgerr(k+1) = (sol - x(:,k+1))'*A*(sol-x(:,k+1))
    end


end