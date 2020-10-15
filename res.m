% function [r,dr] = res(x,v,W,u,fun,dfun,d2fun,d3fun)
% % BVP for the Poisson equation is setup here
% %
% % residual functions r and theor derivatives w.r.t. parameters dr 
% % are evaluated in this function
% %
% % computer diff_operator(Psi(x)) - RHS(x)
% % boundary functions
% [~,~,~,~,d2f0,d2f1,d2g0,d2g1,h,hp,hpp,rhs,~] = setup();
% % differential operator is d^2/dx^2 + d^2/dy^2
% % differential operator applied to A(x,y), the bdry term
% d2A = @(x)(1-x(2))*d2g0(x(1)) + x(2)*d2g1(x(1)) + ...
%     (1-x(1))*d2f0(x(2)) + x(1)*d2f1(x(2));
% % differential operator applied to B(x,y) = x(1-x)y(1-y)NN(x,y,v,W,u)
% h1 = h(x(1));
% h2 = h(x(2));
% hp1 = hp(x(1));
% hp2 = hp(x(2));
% [f,fx,fy,fxx,fyy,df,dfx,dfy,dfxx,dfyy] = NN(x,v,W,u,fun,dfun,d2fun,d3fun);
% d2B = hpp*(h1+h2)*f + 2*(hp1*h2*fx + hp2*h1*fy) + h1*h2*(fxx + fyy);
% % residual r = d2A + d2B - RHS
% r = d2A(x) + d2B - rhs(x(1),x(2));
% % derivative of r w.r.t. parameters
% dr = hpp*(h1+h2)*df + 2*(hp1*h2*dfx + hp2*h1*dfy) + h1*h2*(dfxx + dfyy);
% end

%% mixed boundary condition
function [r,dr] = res(x,v,W,u,fun,dfun,d2fun,d3fun,d4fun)
% BVP for the Poisson equation is setup here
%
% residual functions r and theor derivatives w.r.t. parameters dr 
% are evaluated in this function
%
% computer diff_operator(Psi(x)) - RHS(x)
% boundary functions
[~,~,~,~,d2f0,d2f1,d2g0,d2g1,hx,hxp,hxpp,hy,hyp,hypp,rhs,~] = setup();
% differential operator is d^2/dx^2 + d^2/dy^2
% differential operator applied to A(x,y), the bdry term
d2A = @(x)d2g0(x(1)) + x(2)*d2g1(x(1)) + ...
    (1-x(1))*d2f0(x(2)) + x(1)*d2f1(x(2));
% differential operator applied to B(x,y) = x(1-x)y(1-y)NN(x,y,v,W,u)
h1 = hx(x(1));
h2 = hy(x(2));
hp1 = hxp(x(1));
hp2 = hyp(x(2));
hpp1 = hxpp;
hpp2 = hypp;
[f,fx,fy,fxx,fyy,fyx,fyxx,fyyy,df,dfx,dfy,dfxx,dfyy,dfyx,dfyxx,dfyyy] = NN(x,v,W,u,fun,dfun,d2fun,d3fun,d4fun);
x1 = [x(1);1];
[f1,fx1,fy1,fxx1,fyy1,fyx1,fyxx1,fyyy1,df1,dfx1,dfy1,dfxx1,dfyy1,dfyx1,dfyxx1,dfyyy1] = NN(x1,v,W,u,fun,dfun,d2fun,d3fun,d4fun);

d2B = (h1*hpp2+h2*hpp1)*(f-f1-fy1) + 2*(hp1*h2*(fx-fx1-fyx1) + hp2*h1*(fy-fyy1)) ...
    + h1*h2*(fxx - fxx1 - fyxx1 + fyy - fyyy1);
% residual r = d2A + d2B - RHS
r = d2A(x) + d2B - rhs(x(1),x(2));
% derivative of r w.r.t. parameters
dr = (h1*hpp2+h2*hpp1)*(df-df1-dfy1) + 2*(hp1*h2*(dfx-dfx1-dfyx1) ...
    + hp2*h1*(dfy-dfyy1)) + h1*h2*(dfxx - dfxx1 - dfyxx1 + dfyy - dfyyy1);
end



