function []=plot_error(E,opts)
% Create plots illustrating design optimization and reliability assesment
%
% Inputs:
% E = structure with error model
% opts = structure with inputs and options
%
% Outputs:
% none
%==========================================================================

close all

fontName='Arial';
fontSize=8;

udet=opts.design.udet;

%% Find optimum design
k=0;
[xopt_norm,fopt]=Design_Opt(k,udet,E,opts);
xopt=norm_design(xopt_norm,'x',0,opts);

%% Find MPP
opts.reliability.display='iter'; 
[Pf,B,MPP_norm]=Reliability_Analysis(xopt_norm,E,opts,'low_upd');
MPP=norm_design(MPP_norm,'u',0,opts);
A=[500,1000];
B=[0.2,0.1].*A;
MPP_stand_norm=(MPP-A)./B;

%% Design optimization plot
% Plot variables
n=40;
x1_norm=linspace(0,1,n)';
x2_norm=linspace(0,1,n)';
x=norm_design([x1_norm,x2_norm],'x',0,opts);
x1=x(:,1);
x2=x(:,2);
[X1_grid,X2_grid]=meshgrid(x1_norm,x2_norm);
X1=reshape(X1_grid,n^2,1);
X2=reshape(X2_grid,n^2,1);

% Low-fidelity mean
E.returnKstd=false;
g=g_H([X1,X2],repmat(udet,n^2,1),E,opts);
G.mean=reshape(g,n,n);

% Low-fidelity confidence interval
E.returnKstd=true;
E.k=1.96;
g=g_H([X1,X2],repmat(udet,n^2,1),E,opts);
G.lb=reshape(g,n,n);
E.k=-1.96;
g=g_H([X1,X2],repmat(udet,n^2,1),E,opts);
G.ub=reshape(g,n,n);

% High-fidelity model
g=g_H_true([X1,X2],repmat(udet,n^2,1),opts);
g_true=reshape(g,n,n);

dtip=2.25*(0.1^3)
g_L=dtip-eul_def([X1,X2],repmat(udet,n^2,1),opts);
g_L=reshape(g_L,n,n);
error=E.Ei{1};
G.lb=g_L+min(error);
G.ub=g_L+max(error);

figure()
hold on
[~,h_true]=contour(x1,x2,g_true,[0,0],'-k');
[~,h_mean]=contour(x1,x2,G.mean,[0,0],'-m');
[~,h_lb]=contour(x1,x2,G.lb,[0,0],'--m');
[~,h_ub]=contour(x1,x2,G.ub,[0,0],'--m');
h_opt=plot(xopt(1),xopt(2),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',3);
hh=legend([h_true,h_mean,h_lb,h_opt],'$g_H(\mathbf{x},\mathbf{u}_{det})=0$',...
    '$\bar{g}_H(\mathbf{x},\mathbf{u}_{det})=0$','95\% CI',...
    'Optimum when k=0');
set(hh,'Interpreter','latex')
xlabel('Width (in), x_1')
ylabel('Thickness (in), x_2')

x=[2.75 2.75];
y=[1.75 1.75];
[xa,ya]=ds2nfu(x,y);
dim = [xa(1) ya(1)+0.05 .05 .05];
str='Infeasible';
annotation('textbox',dim,'String',str,'linestyle','none','FontName',fontName,'FontSize',fontSize)

x=[4 4];
y=[3 3];
[xa,ya]=ds2nfu(x,y);
dim = [xa(1) ya(1)+0.05 .05 .05];
str='Feasible';
annotation('textbox',dim,'String',str,'linestyle','none','FontName',fontName,'FontSize',fontSize)

%% Design optimization plot
% Plot variables
n=40;
u1_stand_norm=linspace(0,6,n)';
u2_stand_norm=linspace(0,7,n)';
U_norm=myTinv([u1_stand_norm,u2_stand_norm],opts);

u1_norm=U_norm(:,1);
u2_norm=U_norm(:,1);
u=norm_design([u1_norm,u2_norm],'u',0,opts);
u1=u(:,1);
u2=u(:,2);
[U1_grid,U2_grid]=meshgrid(u1_norm,u2_norm);
U1=reshape(U1_grid,n^2,1);
U2=reshape(U2_grid,n^2,1);

% Low-fidelity mean
E.returnKstd=false;
g=g_H(repmat(xopt_norm,n^2,1),[U1,U2],E,opts);
G.mean=reshape(g,n,n);

% Low-fidelity confidence interval
E.returnKstd=true;
E.k=1.96;
g=g_H(repmat(xopt_norm,n^2,1),[U1,U2],E,opts);
G.lb=reshape(g,n,n);
E.k=-1.96;
g=g_H(repmat(xopt_norm,n^2,1),[U1,U2],E,opts);
G.ub=reshape(g,n,n);

% High-fidelity model
g=g_H_true(repmat(xopt_norm,n^2,1),[U1,U2],opts);
g_true=reshape(g,n,n);

dtip=2.25*(0.1^3)
g_L=dtip-eul_def(repmat(xopt_norm,n^2,1),[U1,U2],opts);
g_L=reshape(g_L,n,n);
error=E.Ei{1};
G.lb=g_L+min(error);
G.ub=g_L+max(error);

figure()
hold on
[~,h_true]=contour(u1_stand_norm,u2_stand_norm,g_true,[0,0],'-k');
[~,h_mean]=contour(u1_stand_norm,u2_stand_norm,G.mean,[0,0],'-m');
[~,h_lb]=contour(u1_stand_norm,u2_stand_norm,G.lb,[0,0],'--m');
[~,h_ub]=contour(u1_stand_norm,u2_stand_norm,G.ub,[0,0],'--m');
% h_opt=plot(MPP_stand_norm(1),MPP_stand_norm(2),'o','MarkerEdgeColor','k','MarkerFaceColor','r');
hh=legend([h_true,h_mean,h_lb],'$g_H(\mathbf{x}_{opt},\mathbf{U})=0$',...
    '$\bar{g}_H(\mathbf{x}_{opt},\mathbf{U})=0$','95\% CI');
set(hh,'Interpreter','latex')

x=[0.5 0.5];
y=[0.5 0.5];
[xa,ya]=ds2nfu(x,y);
dim = [xa(1) ya(1)+0.05 .05 .05];
str='Safe';
annotation('textbox',dim,'String',str,'linestyle','none','FontName',fontName,'FontSize',fontSize)

x=[4 4];
y=[4 4];
[xa,ya]=ds2nfu(x,y);
dim = [xa(1) ya(1)+0.05 .05 .05];
str='Failure';
annotation('textbox',dim,'String',str,'linestyle','none','FontName',fontName,'FontSize',fontSize)

xlabel('Horizontal load ($\sigma$), $\hat{u}_1$','Interpreter','latex')
ylabel('Vertical load ($\sigma$), $\hat{u}_1$','Interpreter','latex')

end