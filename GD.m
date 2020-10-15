function [fall,norg] = GD(nt,N,tol,iter_max)
fsz = 16; % fontsize
%%
eta = 0.1;
gam = 0.9;
% iter_max = 100;
% tol = 5e-3;
jmax = ceil(log(1e-14)/log(gam)); % max # of iterations in line search
%% setup training mesh
% nt = 5;
t = linspace(0,1,nt+2);
[xm,ym] = meshgrid(t,t);
I = 2:(nt+1);
xaux = xm(I,I);
yaux = ym(I,I);
xy = [xaux(:),yaux(:)]';
%% initial guess for parameters
% N = 10; % the number of hidden nodes
npar = 4*N;
w = ones(npar,1);
%%
[r,J] = Res_and_Jac(w,xy);
f = F(r);
g = J'*r;
nor = norm(g);

fprintf('Initially: f = %d, nor(g) = %d\n',f,nor); 
%% The trust region BFGS method
tic

iter = 0;
I = eye(length(w));
% quadratic model: m(p) = (1/2)||r||^2 + p'*J'*r + (1/2)*p'*J'*J*p;
norg = zeros(iter_max+1,0);
fall = zeros(iter_max+1,0);
norg(1) = nor;
fall(1) = f;
while nor > tol && iter < iter_max
    % solve the constrained minimization problem using dogleg strategy
    p = -g;
    a = 1;
    aux = eta*g'*p;
    for j = 0 : jmax
        wtry = w + a*p;
        [rtry, Jtry] = Res_and_Jac(wtry,xy);
        f1 = F(rtry);
        if f1 < f + a*aux
            break;
        else
            a = a*gam;
        end
    end
    w = w + a*p;
    [r,J] = Res_and_Jac(w,xy);
    f = F(r);
    g = J'*r;
    nor = norm(g);
    fprintf('iter %d: line search: j = %d, a = %d, f = %d, norg = %d\n',iter,j,a,f,nor);
    iter = iter + 1;
    norg(iter+1) = nor;
    fall(iter+1) = f;
end
fprintf('iter # %d: f = %.14f, |df| = %.4e\n',iter,f,nor);
cputime = toc;
fprintf('CPUtime = %d, iter = %d\n',cputime,iter);
%% visualize the solution
nt = 101;
t = linspace(0,1,nt);
[xm,ym] = meshgrid(t,t);
[fun,dfun,~,~] = ActivationFun();
[v,W,u] = param(w);
[f0,f1,g0,g1,~,~,~,~,hx,~,~,hy,~,~,~,exact_sol] = setup();
% A = @(x,y)(1-x).*f0(y) + x.*f1(y) + (1-y).*(g0(x)-((1-x)*f0(0)+x*f1(0))) + ...
%      y.*(g1(x)-((1-x)*f0(1)+x*f1(1)));
% B = h(xm).*h(ym);
A = @(x,y)(1-x).*f0(y) + x.*f1(y) + g0(x) - ((1-x).*g0(0)+x.*g0(1)) ...
    + y.*(g1(x)-((1-x).*g1(0) + x.*g1(1)));
B = hx(xm).*hy(ym);
NNfun = zeros(nt);
for i = 1 : nt
    for j = 1 : nt
        x = [xm(i,j);ym(i,j)];
        NNfun(i,j) = v'*fun(W*x + u);
    end
end
NNfun1 = zeros(nt);
for i = 1 : nt
    for j = 1 : nt
        x = [xm(i,j);1];
        NNfun1(i,j) = v'*fun(W*x + u);
    end
end
NNfun1dy = zeros(nt);
for i = 1 : nt
    for j = 1 : nt
        x = [xm(i,j);1];
        NNfun1dy(i,j) = v'*(W(:,2).*dfun(W*x + u));
    end
end
sol = A(xm,ym) + B.*(NNfun-NNfun1-NNfun1dy);
esol = exact_sol(xm,ym);
err = sol - esol;
fprintf('max|err| = %d, L2 err = %d\n',max(max(abs(err))),norm(err(:)));

%
folder = 'figs_gd/';
figure(1);clf;
contourf(t,t,sol,linspace(min(min(sol)),max(max(sol)),20));
colorbar;
set(gca,'Fontsize',fsz);
xlabel('x','Fontsize',fsz);
ylabel('y','Fontsize',fsz);
filename = [folder,'sol','.png'];
saveas(gcf,filename)
%
figure(2);clf;
contourf(t,t,err,linspace(min(min(err)),max(max(err)),20));
colorbar;
set(gca,'Fontsize',fsz);
xlabel('x','Fontsize',fsz);
ylabel('y','Fontsize',fsz);
filename = [folder,'err','.png'];
saveas(gcf,filename)
%
figure(3);clf;
subplot(2,1,1);
fall(iter+2:end) = [];
plot((1:iter+1)',fall,'Linewidth',2,'Marker','.','Markersize',20);
grid;
set(gca,'YScale','log','Fontsize',fsz);
xlabel('k','Fontsize',fsz);
ylabel('f','Fontsize',fsz);
subplot(2,1,2);
norg(iter+2:end) = [];
plot((1:iter+1)',norg,'Linewidth',2,'Marker','.','Markersize',20);
grid;
set(gca,'YScale','log','Fontsize',fsz);
xlabel('k','Fontsize',fsz);
ylabel('|| grad f||','Fontsize',fsz);
filename = [folder,'f_grad','.png'];
saveas(gcf,filename)
end

%%
function p = cauchy_point(B,g,R)
    ng = norm(g);
    ps = -g*R/ng;
    aux = g'*B*g;
    if aux <= 0
        p = ps;
    else
        p = min(ng^3/(R*aux),1);
    end
end
%%
function f = F(r)
    f = 0.5*r'*r;
end
