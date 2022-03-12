function [] = smootheLB(F,V,nharm,f,savenamef,savenameB)
% f denotes function values to smoothe
% byuIn gives surface that LB will get

addpath /cis/home/kstouff4/Functions/byuFunctionsDT
addpath /cis/home/kstouff4/Documents/SurfaceTools/FromDaniel/younesShapeFun/Surfaces

if (min(F) < 1)
    F = F+1;
end
[B,D,L,A] = computeHarmonics(F,V,nharm,1); % boundary condition 1 is von Neuman like in 2005 Anqi paper, B is eigenvectors, D is eigenvalues
disp(max(max(L*B - A*B*D)))
disp(min(min(L*B - A*B*D)))
disp(B(1:4,1:4))
disp(D(1:4,1:4))
A = full(A);
L = full(L);
Ah = A^(1/2);
Bh = Ah*B;
disp(max(max(L*Bh - Bh*A*D)))
disp(min(min(L*Bh - Bh*A*D)))
disp(max(max(transpose(Bh)*Bh)))
disp(min(min(transpose(Bh)*Bh)))
save(savenameB,'B','D','A','L','Ah','Bh');

if (f ~= 0)
    newf = transpose(B)*f; % should be nharm x 1 i.e. gives coefficients in new basis degree x numV * numV x numF = degree x numF
    newf = B*newf; % represent function as sum of coeff*basis vectors numV x degree * degree x numF
    newfL = inv(L)*f;
    save(savenamef,'newf','newfL');
end
