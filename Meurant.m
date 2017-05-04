%This script function plots a Visualization of the Theorem by Meurant 1999
%which relates the Gauss Christoffel quadrature to the error in the CG
%iteration
%
% Input: 
% A : SPD nxn matrix 
% b : nx1 vector 
% increments: iterations at which the plot will be shown

function [ ] = Meurant(A,b,increments)

    n = size(A,1);
    [U,S] = eig(full(A));
    S = diag(S);

    [~, cgerr] = VanillaCG(A,b);
    cgerr = cgerr/(norm(b)^2);

    v1 = b/norm(b);

    
    y = (v1'*U).^2;
    y = y(:);

    sy = cumsum(y);
    omega = sy;

    lambda = (S(1)):0.001:(S(end)*1.1);
    omegafun = @(lambda) ((S < lambda).*(lambda <= [S(2:end); inf]))'*omega;

    for i=1:length(lambda)
        om(i) = omegafun(lambda(i));
    end

    
    
    clf;
    Id = lambda;
    f = ones(size(lambda,1))./lambda;


    % Plot of approximation of omega
    subplot(1,3,1);
    
    hold on;
    X = [om,  zeros(size(om))];
    Y = [Id, fliplr(Id)];
    %area(Id,om,'Facecolor','y');
    plot(Id,om,'b','Linewidth',2);
    plot(S,zeros(size(omega)),'bo');
    
    
    xlim([0,max(S)]);
    ylim([0,max(om)]);
    ylabel('$$ \omega(\lambda) $$','Interpreter','Latex','Fontsize',20);
    xlabel('$$\lambda$$','Interpreter','Latex','Fontsize',20);
    set(gca,'fontsize',10)



    subplot(1,3,2);

    X = [om,  zeros(size(om))];
    Y = [f, fliplr(f)];
    hold on;
    area(om,f,'FaceColor','y');

    %fill(X,Y,'c'):

    


    xlim([0,1]);
    ylim([0,1.1*max(f.*(om>0))]);
    plot(om,f,'b','Linewidth',2);
    plot(omega,zeros(size(omega)),'bo');
    xlabel('$$\omega(\lambda)$$','Interpreter','Latex','Fontsize',20);
    ylabel({'$$f(\lambda)= \frac{1}{\lambda}$$'},'Interpreter','Latex','Fontsize',20);
    set(gca,'fontsize',10)
    
    pbaspect([1 max(f.*(om>0)) 1])


    %Plot of CG error
    subplot(1,3,3);
    semilogy(0,cgerr(1),'ro', 'Markerfacecolor', 'r' );
    xlabel('Iteration','Interpreter','Latex','Fontsize',20);
    ylabel('$$ \frac{\| x - x_k \|_A^2}{\|r_0\|_2^2}$$','Interpreter','Latex','Fontsize',20);
    set(gca,'yscale','log');
    set(gca,'fontsize',10)
    ylim([min(cgerr), max(cgerr)]);
    xlim([0,1]);
    
    
    pause()



    % Start Lanczos: Stolen from Dr. Bindel's lecture notes for CS6220
    [Q,alpha, beta]=lanczos(A,b,n);
    T = diag(alpha) + diag(beta(1:end-1),1) + diag(beta(1:end-1),-1);





    for m =1:length(increments)
        k = increments(m);
        [evecs,Sk] = eig(T(1:k,1:k));
        Sk = diag(Sk);
        jumps = evecs(1,:)'.^2;
        omegahat = cumsum(jumps);
        omegahatfun= @(lambda) ((Sk < lambda).*(lambda <= [Sk(2:end); inf]))'*omegahat;

        for i=1:length(lambda)
            omhat(i) = omegahatfun(lambda(i));
        end
        %close all;

        clf

        %fill(X,Y,'y');   
        [pos,neg] = get_polygons(om, omhat);

        % Plot of approximation of omega
        subplot(1,3,1);
        %clf
        hold on;


        plot(Id,om,'b','Linewidth',2);
        plot(Id,omhat,'r','Linewidth',2);
        plot(S,zeros(size(S)),'bo');
        plot(Sk,zeros(size(Sk)),'ro');
        
       % legend('$$ \omega $$','$$ \hat{\omega}(\lambda) $$','Interpreter','Latex')
        xlim([0,max(S)]);
        ylim([0,max(om)]);

        ylabel('$$ \omega(\lambda) $$','Interpreter','Latex','Fontsize',20);
        xlabel('$$\lambda$$','Interpreter','Latex','Fontsize',20);
        set(gca,'fontsize',10)




        subplot(1,3,2);
        %fill(X,Y,'y');   
        [pos,neg] = get_polygons(om, omhat);

        hold on;
        for i =1:length(neg)
           X = [om(neg{i}), fliplr(omhat(neg{i}))];
           Y = [f(neg{i}), fliplr(f(neg{i}))];
           fill(X,Y,'y');
        end

        for i =1:length(pos)
           X = [om(pos{i}), fliplr(omhat(pos{i}))];
           Y = [f(pos{i}), fliplr(f(pos{i}))];
           fill(X,Y,'c');
        end



        xlim([0,1]);
        ylim([0,1.1*max(f.*(om>0))]);
        plot(om,f,'b','Linewidth',2);
        plot(omhat,f,'r','Linewidth',2);
        plot(omega,zeros(size(omega)),'bo');
        plot(omegahat,zeros(size(omegahat)),'ro');
        xlabel('$$\omega(\lambda)$$','Interpreter','Latex','Fontsize',20);
        ylabel({'$$f(\lambda)= \frac{1}{\lambda}$$'},'Interpreter','Latex','Fontsize',20);
       % legend({'[$$ \omega $$','$$ \hat{\omega}(\lambda) $$'},'Interpreter','Latex')
        set(gca,'fontsize',10)

        pbaspect([1 max(f.*(om>0)) 1])

        
        %Plot of CG error
        subplot(1,3,3);
        hold on;
        semilogy(0:k-1, cgerr(1:k),'ro')
        
        semilogy(k,cgerr(k+1),'ro', 'Markerfacecolor', 'r' );
        xlabel('Iteration','Interpreter','Latex','Fontsize',20);
        ylabel('$$ \frac{\| x - x_k \|_A^2}{\|r_0\|_2^2}$$','Interpreter','Latex','Fontsize',20);
        set(gca,'yscale','log');
        set(gca,'fontsize',10)

        ylim([min(cgerr), max(cgerr)]);
        xlim([0,k+1]);


        %axis equal

        shg
        pause()
    
    
end