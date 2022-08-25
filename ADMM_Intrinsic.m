

function [z,x,y, error] = ADMM_Intrinsic(c,niter,rho, lmb)

n = size(c,3);
mx = size(c,2);
my = size(c,1);

vc = zeros(mx*my,n);
for i = 1:n
    tmp = c(:,:,i);
    vc(:,i) = tmp(:);
end

S1 = sparse(1:mx*my,1:mx*my,-2*ones(1,mx*my),mx*my,mx*my);
S2 = sparse(2:mx*my,1:mx*my-1,ones(1,mx*my-1),mx*my,mx*my);
S = S2 + S1 + S2';
StS = S'*S;
StSrho = lmb*StS + rho*speye(mx*my);

e = ones(mx*my,1);
Dx = spdiags([e -e], 0:1, mx*my,mx*my);
DxtDx = Dx'*Dx;
Dy = Dx';
DytDy = Dy'*Dy;
DtD = DxtDx + DytDy;
DtDrho = DtD + rho*speye(mx*my);
DtDrhovc = DtDrho*vc;

z = -ones(mx*my,1);
x = -ones(mx*my,n);
y = -ones(mx*my,n);

error = nan(1,niter);
for k = 1:niter
    
    x = DtDrho \ ( DtDrhovc - rho*z - y );
    x(:,1) = x(:,1).*double(x(:,1) <= 0);
    x(:,2) = x(:,2).*double(x(:,2) <= 0);
    x(:,3) = x(:,3).*double(x(:,3) <= 0);
    
    z = (n*StSrho)\sum(rho*vc - rho*x - y,2);
    z = z.*double(z <= 0);
    
    y = y + rho*(x + z - vc);
        
end


