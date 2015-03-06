%% Function used only to plot the error

function [error] = errorfunction(Nr)
    %% Mesh
    [p, tri, edge] = getDisk(Nr);

    %% Assembly
    [Ah] = stiffness_2D(Nr,p,tri);
    [fh] = load_2D(Nr,p,tri);

    %% Imposing the dirichlet boundary
    index = length(edge(:,1)) - 1;
    for i=(Nr-index):Nr
        fh(i)=0;
        Ah(i,:)=0;
        Ah(:,i)=0;
        Ah(i,i)=1;
    end

    %% Numerical Solution
    uh = Ah\fh;

    %% Analytical Solution
    ua=@(x,y) sin(2*pi*(x^2+y^2));
    uanal=zeros(Nr,1);
    for i=1:Nr
        uanal(i)=ua(p(i,1),p(i,2));
    end

    %% Error
    error = max(abs(uh-uanal));
end