
clear
close all


original = double(imread('turtle.png'));

obj(:,:,1) = original(:,:,1)/max(max(original(:,:,1)));
obj(:,:,2) = original(:,:,2)/max(max(original(:,:,2)));
obj(:,:,3) = original(:,:,3)/max(max(original(:,:,3)));

c = zeros(size(obj));
for k = 1:size(obj,3)
    tmp = log(obj(:,:,k));
    indx = find( tmp == -inf);
    tmp(indx) = 0;
    c(:,:,k) = tmp;
end

tic
niter = 1500;
rho = 1e-2;
lmb = 1e2;
[z,x,y, error] = ADMM_Intrinsic(c,niter,rho, lmb);

S = zeros(size(z));
indx = find(z(:));
S(indx) = exp(z(indx));
S = reshape(S, size(c,1), size(c,2));

figure, imagesc(S), axis image, colormap gray, title('Shape')

indx1 = find(S(:) == 0);
S(indx1) = 1;
figure, imagesc(0.6*obj./S), axis image, title('Reflectance')

