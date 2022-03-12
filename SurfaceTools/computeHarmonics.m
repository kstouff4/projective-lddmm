function [B, D, L, A] = computeHarmonics(F,V, nharm, bCond)
%%% Data processing
%%
%% Computes nharm eigenvectors with smallest eigenvalues in absolute value
%% for the Laplacian operator
%% fv is the triangulated surface
%% nharm is the number of eigenvectors to compute
%% bCond specifies the boundary condition: default is Dirichlet; use 1 for von Neuman.
%% [B,D] are eigenvectors and eigenvalues as returned by eigs

if (nargin == 2)
    bCond = 1;
end ;

[L, A, bVert] = laplacianMatrix(F,V) ;

I = find(bVert == 0) ;
J = find(bVert==1) ;
if (isempty(J)) %closed surface
    [B,D] = eigs(L, A, size(A,1), 1) ;
    %D = diag(D) ;
    disp('closed surface')
else if (bCond == 1) %von Neuman boundary condition 
        %opts.tol=1e-10;
        LI = L  ;
        AI = sparse(size(L,1), size(L,2)) ;
        for k=1:numel(I),
            AI(I(k), I(k)) = A(I(k), I(k)) ;
        end
        [B, D] = eigs(LI, AI, nharm, 1) ;
        D = diag(D) ;
    else % Dirichlet boundary condition
        L0 = L(I, I) ;
        A0 = A(I, I) ;
        [B0,D] = eigs(L0, A0, nharm, 1) ;
        D = diag(D) ;
        B = zeros(size(L,1), nharm) ;
        B(I, :) = B0 ;
    end
end
