function [L, A, bVert] = laplacianMatrix(F,V)
%%% [L, A,bVert] = laplacianMatrix(fv)
%% Computes Laplacian operator on triangulated surface fv
%% Returns L as a sparse matrix
%% A is a diagonal sparse matrix that contains areas attributed to vertices
%% bVert is a flag on vertices (1 if boundary, 0 otherwise)


%F = fv.faces ;
%V = fv.vertices ;
nf = size(F,1) ;
nv = size(V,1) ;
%% compute areas of faces and vertices
AF = zeros(nf,1) ;
AV = zeros(nv, 1) ;
for k=1:nf,
    %% determining if face is obtuse
    x12 = V(F(k,2), :) - V(F(k,1), :) ;
    x13 = V(F(k,3), :) - V(F(k,1), :) ;
    n12 = norm(x12) ;
    n13 = norm(x13) ;
    c1 = sum(sum(x12.*x13))/(n12*n13) ;
    x23 = V(F(k,3), :) - V(F(k,2), :) ;
    n23 = norm(x23) ;
    c2 = - sum(sum(x12 .*x23)) / (n12 * n23) ;
    c3 = sum(sum(x13.*x23)) / (n13*n23) ;
    AF(k) = norm(cross(x12, x13))/2 ;
    if (c1 < 0)  %face obtuse at vertex 1
        AV(F(k,1)) = AV(F(k,1)) + AF(k)/2 ;
        AV(F(k,2)) = AV(F(k,2)) + AF(k)/4 ;
        AV(F(k,3)) = AV(F(k,3)) + AF(k)/4 ;
    else if (c2 < 0) %face obuse at vertex 2
        AV(F(k,1)) = AV(F(k,1)) + AF(k)/4 ;
        AV(F(k,2)) = AV(F(k,2)) + AF(k)/2 ;
        AV(F(k,3)) = AV(F(k,3)) + AF(k)/4 ;
        else if (c3 < 0) % face obtuse at vertex 3
        AV(F(k,1)) = AV(F(k,1)) + AF(k)/4 ;
        AV(F(k,2)) = AV(F(k,2)) + AF(k)/4 ;
        AV(F(k,3)) = AV(F(k,3)) + AF(k)/2 ;
            else %non obtuse face
                cot1 = c1 / sqrt(1-c1^2) ;
                cot2 = c2 / sqrt(1-c2^2) ;
                cot3 = c3 / sqrt(1-c3^2) ;
                AV(F(k,1)) = AV(F(k,1)) + (sum(sum(x12.^2)) * cot3 + sum(sum(x13.^2)) * cot2)/8 ;
                AV(F(k,2)) = AV(F(k,2)) + (sum(sum(x12.^2)) * cot3 + sum(sum(x23.^2)) * cot1)/8 ;
                AV(F(k,3)) = AV(F(k,3)) + (sum(sum(x13.^2)) * cot2 + sum(sum(x23.^2)) * cot1)/8 ;
            end
        end
    end
end

for k=1:nv,
    if (abs(AV(k)) <1e-10)
        disp(sprintf('Warning: vertex %d has no face; use removeIsolated.m', k)) ;
    end
end

%% compute edges and detect boundary
edm = sparse(nv,nv) ;
E = zeros(3*nf, 2) ;
j = 0 ;
for k=1:nf,
    if (edm(F(k,1), F(k,2))== 0)
        j = j+1 ;
        edm(F(k,1), F(k,2)) = j ;
        edm(F(k,2), F(k,1)) = j ;
        E(j, :) = [F(k,1), F(k,2)] ;
    end
    if (edm(F(k,2), F(k,3))== 0)
        j = j+1 ;
        edm(F(k,2), F(k,3)) = j ;
        edm(F(k,3), F(k,2)) = j ;
        E(j, :) = [F(k,2), F(k,3)] ;
    end
    if (edm(F(k,1), F(k,3))== 0)
        j = j+1 ;
        edm(F(k,3), F(k,1)) = j ;
        edm(F(k,1), F(k,3)) = j ;
        E(j, :) = [F(k,3), F(k,1)] ;
    end
end
E = E(1:j, :) ;
edgeFace = sparse(j, nf) ;
ne = j ;
for k=1:nf,
    edgeFace(edm(F(k,1), F(k,2)), k) = 1 ;
    edgeFace(edm(F(k,2), F(k,3)), k) = 1 ;
    edgeFace(edm(F(k,3), F(k,1)), k) = 1 ;
end
    
bEdge = zeros(ne, 1) ;
bVert = zeros(nv, 1) ;
edgeAngles = zeros(ne, 2) ;
for k=1:ne,
    I = find(edgeFace(k, :) == 1) ;
    for u=1:numel(I),
        f = I(u) ;
        i1 = find(F(f, :) == E(k,1)) ;
        i2 = find(F(f, :) == E(k,2)) ;
        s = i1+i2 ;
        if (s == 3)
            i3 = 3 ;
        else if (s==4)
                i3 = 2 ;
            else if (s==5)
                    i3=1 ;
                end
            end
        end
        x1 = V(F(f,i1), :) - V(F(f,i3), :) ;
        x2 = V(F(f,i2), :) - V(F(f,i3), :) ;
        a = sum(sum(cross(x1, x2) .* cross(V(F(f,2), :) - V(F(f,1), :), V(F(f, 3), :) - V(F(f, 1), :)))) ;
        b = sum(sum(x1.*x2)) ;
        if (a  > 0)
            edgeAngles(k, u) = b/sqrt(a) ;
        else
            edgeAngles(k, u) =  b/sqrt(-a) ;
        end
    end
    if (numel(I) == 1)
        %% boundary edge
        bEdge(k) = 1 ;
        bVert(E(k,1)) = 1 ;
        bVert(E(k,2)) = 1;
        edgeAngles(k,2) = 0; 
    end
end
        

%% Compute Laplacian matrix
L = sparse(nv, nv) ;

for k=1:ne,
    L(E(k,1), E(k,2)) = (edgeAngles(k,1) + edgeAngles(k,2)) /2 ;
    L(E(k,2), E(k,1)) = L(E(k,1), E(k,2)) ;
end

for k=1:nv,
    L(k,k) = - sum(sum(L(k, :))) ;
end

A = sparse(nv, nv) ;
for k=1:nv,
     A(k, k) = AV(k) ;
end ;
        