function [] = PlotHists(dat,opts)
% Generate histograms from MCS data

%% Unpack variables
v2struct(dat);
dim_x=opts.design.dim;
dim_u=opts.reliability.dim;

%% Define plot options
fc1='b';            % Color of initial histograms
fc2='y';            % Color of final histograms
ec='k';             % Edge color of histograms
falph=0.4;          % Transparency of histograms
nbins=60;           % Number of bins in histograms
ms=3;               % Marker size

% Initialize figure counter
numFig=1;           

% Set graphics renderer
opengl software

%% Design variable distribution
if ~opts.fig.saveFiles
    figure(numFig)
    numFig=numFig+1;
    clf
end
d=opts.design.dim;
r=2;
c=ceil(d/r);
lab=opts.fig.xlabels;
for k=1:d
    if d>1 && ~opts.fig.saveFiles
        subplot(r,c,k)
    else
        figure(numFig)
        numFig=numFig+1;
        clf
    end
    hold on

    if isfield(opts.fig,'xlims')
        x_edge=linspace(opts.fig.xlims.x(k,1),opts.fig.xlims.x(k,2),nbins+1);
    else
        x_edge=linspace(norm_design(min(min(xini(:,k),Xfinal(:,k))),...
            num2str(k),0,opts),norm_design(max(max(xini(:,k),Xfinal(:,k))),num2str(k),0,opts),nbins+1);
    end        

    x_edge=[x_edge(1:nbins),inf];

    y=histc(norm_design(Xfinal(:,k),num2str(k),0,opts),x_edge);
    h2=bar(x_edge(1:nbins),y(1:nbins),'histc');
    set(h2,'facecolor',fc2,'EdgeColor',ec,'linestyle','-')

    h11=plot(norm_design(xini(:,k),num2str(k),0,opts),0,'s','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','b');
    h22=plot(mean(norm_design(Xfinal(:,k),num2str(k),0,opts)),0,'o','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','g');

    [~,BLicons]=legend([h2 h11 h22],'Final','Initial Design','Final Mean');
    if isfield(opts.fig,'labels')
        xlabel(opts.fig.labels{1})
    else
        xlabel([lab{k},', x_',num2str(k)])
    end
    ylabel('Number of Occurrences')

    set(findobj(gca,'type','patch'),'facea',falph)
    set(findobj(BLicons, 'type', 'patch'),'facea',falph)

    if opts.fig.saveFiles
        fig_export(['x',num2str(k),opts.fig.fname],0.33)
        [~,BLicons]=legend('-DynamicLegend');
        set(findobj(BLicons, 'type', 'patch'),'facea',falph)
        fig_export(['x',num2str(k),opts.fig.fname],0.33)
    end
end

%% Objective distribution
figure(numFig)
numFig=numFig+1;
clf
hold on

if isfield(opts.fig,'xlims')
    x_edge=linspace(opts.fig.xlims.f(1,1),opts.fig.xlims.f(1,2),nbins+1);
else
    x_edge=linspace(min(min(Fini,Ffinal)),max(max(Fini,Ffinal)),nbins+1);
end   
x_edge=[x_edge(1:nbins),inf];

y=histc(Fini,x_edge);
h1=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h1,'facecolor',fc1,'EdgeColor',ec,'linestyle','-')

y=histc(Ffinal,x_edge);
h2=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h2,'facecolor',fc2,'EdgeColor',ec,'linestyle','-')

h11=plot(mean(Fini),0,'s','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','b');
h22=plot(mean(Ffinal),0,'o','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','g');

[~,BLicons]=legend([h1 h2 h11 h22],'Initial','Final','Initial Mean','Final Mean');
if isfield(opts.fig,'flabel')
    xlabel([opts.fig.flabel])
else
    xlabel('f(x)')
end
ylabel('Number of Occurrences') 

set(findobj(gca,'type','patch'),'facea',falph)
set(findobj(BLicons, 'type', 'patch'),'facea',falph)

if isfield(opts.fig,'ylims')
    ylim([opts.fig.ylims.f(1,1),opts.fig.ylims.f(1,2)])
end

if opts.fig.saveFiles
    fig_export(['obj',opts.fig.fname])
    [~,BLicons]=legend('-DynamicLegend');
    set(findobj(BLicons, 'type', 'patch'),'facea',falph)
    fig_export(['obj',opts.fig.fname])
end

%% Safety margin distribution
figure(numFig)
numFig=numFig+1;
clf
hold on

if isfield(opts.fig,'xlims')
    x_edge=linspace(opts.fig.xlims.margins(1,1),opts.fig.xlims.margins(1,2),nbins+1);
else
    x_edge=linspace(min(min(Nmeas_1,Nfinal)),max(max(Nmeas_1,Nfinal)),nbins+1);
end   
x_edge=[x_edge(1:nbins),inf];

y=histc(Nmeas_1,x_edge);
h1=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h1,'facecolor',fc1,'EdgeColor',ec,'linestyle','-')

y=histc(Nfinal,x_edge);
h1=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h1,'facecolor',fc2,'EdgeColor',ec,'linestyle','-')

% Redesign rules
xl=xlim;
yl=ylim;
h11=plot(nini,0,'s','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','b');
h22=plot(mean(nre(dat.Q==1)),0,'o','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','c');

if isfield(opts.fig,'ylims')
    ylim([opts.fig.ylims.margins(1,1),opts.fig.ylims.margins(1,2)])
end

plot([nlb nlb],ylim,'--k')
plot([nub nub],ylim,'--k')

plot([0,0],ylim,'--r');

% Legend, labels, and export
[~,BLicons]=legend('Initial','Final','n_{ini}','Mean n_{re}','n_{lb}','n_{ub}','n=0');
xlabel('Safety margin')
ylabel('Number of Occurrences')

set(findobj(gca,'type','patch'),'facea',falph)
set(findobj(BLicons, 'type', 'patch'),'facea',falph)

if opts.fig.saveFiles
    fig_export(['margins',opts.fig.fname])
    [~,BLicons]=legend('-DynamicLegend');
    set(findobj(BLicons, 'type', 'patch'),'facea',falph)
    fig_export(['margins',opts.fig.fname])
end

%% Probability of failure distribution
figure(numFig)
numFig=numFig+1;
clf
hold on
x_edge=linspace(min(min(Pfini,Pffinal)),max(max(Pfini,Pffinal)),nbins+1);
%     x_edge=linspace(0,1,nbins+1);

y=histc(Pfini,x_edge);
h1=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h1,'facecolor',fc1,'EdgeColor',ec,'linestyle','-')

y=histc(Pffinal,x_edge);
h2=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h2,'facecolor',fc2,'EdgeColor',ec,'linestyle','-')

h3=plot(mean(Pfini),0,'s','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','b');
h4=plot(mean(Pffinal),0,'o','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','g');

[~,BLicons]=legend([h1,h2,h3,h4],'Initial','Final','Initial Mean','Final Mean');
xlabel('Probability of Failure')
ylabel('Number of Occurrences') 

set(findobj(gca,'type','patch'),'facea',falph)
set(findobj(BLicons, 'type', 'patch'),'facea',falph)

if opts.fig.saveFiles
    fig_export(['pf',opts.fig.fname])
    [~,BLicons]=legend('-DynamicLegend');
    set(findobj(BLicons, 'type', 'patch'),'facea',falph)
    fig_export(['pf',opts.fig.fname])
end

%% Reliability index distribution
figure(numFig)
numFig=numFig+1;
clf
hold on
if isfield(opts.fig,'xlims')
    x_edge=linspace(opts.fig.xlims.beta(1,1),opts.fig.xlims.beta(1,2),nbins+1);
else
    x_edge=linspace(min(min(Bini,Bfinal)),max(max(Bini,Bfinal)),nbins+1);
end   
x_edge=[x_edge(1:nbins),inf];


y=histc(Bini,x_edge);
h1=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h1,'facecolor',fc1,'EdgeColor',ec,'linestyle','-')

y=histc(Bfinal,x_edge);
h2=bar(x_edge(1:nbins),y(1:nbins),'histc');
set(h2,'facecolor',fc2,'EdgeColor',ec,'linestyle','-')

h3=plot(mean(Bini),0,'s','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','b');
h4=plot(mean(Bfinal),0,'o','MarkerSize',ms,'MarkerEdgeColor','k','MarkerFaceColor','g');

if isfield(opts.fig,'ylims')
    ylim([opts.fig.ylims.beta(1,1),opts.fig.ylims.beta(1,2)])
end

h5=plot([-norminv(opts.margins.pf_bar),-norminv(opts.margins.pf_bar)],ylim,'--r');

[hh,BLicons]=legend([h1,h2,h3,h4,h5],'Initial','Final','Initial Mean','Final Mean','\beta^{\ast}');
set(hh,'Interpreter','latex')
xlabel('Reliability Index')
ylabel('Number of Occurrences') 

set(findobj(gca,'type','patch'),'facea',falph)
set(findobj(BLicons, 'type', 'patch'),'facea',falph)

if opts.fig.saveFiles
    fig_export(['beta',opts.fig.fname])
    [~,BLicons]=legend('-DynamicLegend');
    set(findobj(BLicons, 'type', 'patch'),'facea',falph)
    fig_export(['beta',opts.fig.fname])
end

if d==2
    %% Design variable scatter histogram
    figure(numFig)
    numFig=numFig+1;
    clf
    hold on
    scatterhist(norm_design(Xfinal(:,1),'1',0,opts),norm_design(Xfinal(:,2),'2',0,opts),20);
    xlabel('x_1')
    ylabel('x_2')
    
    %% Design variable 3-D histogram
    figure(numFig)
    numFig=numFig+1;
    clf
    hold on
    if isfield(opts.fig,'xlims')
        x_edge1=linspace(opts.fig.xlims.x(1,1),opts.fig.xlims.x(1,2),nbins+1);
        x_edge2=linspace(opts.fig.xlims.x(2,1),opts.fig.xlims.x(2,2),nbins+1);
    else
        x_edge1=linspace(norm_design(min(min(xini(:,1),Xfinal(:,1))),num2str(1),0),norm_design(max(max(xini(:,1),Xfinal(:,1))),num2str(1),0),nbins+1);
        x_edge2=linspace(norm_design(min(min(xini(:,2),Xfinal(:,2))),num2str(2),0),norm_design(max(max(xini(:,2),Xfinal(:,2))),num2str(2),0),nbins+1);
    end        

    x_edge1=[x_edge1(1:nbins),inf];
    x_edge2=[x_edge2(1:nbins),inf];

    edges{1}=x_edge1;
    edges{2}=x_edge2;
    hist3(norm_design(Xfinal,'x',0,opts),'Edges',edges)
    h = get(gca,'child');
    set(h,'FaceColor','interp', 'CdataMode', 'auto');
    heights=get(h,'Zdata');
    mask=~logical(filter2(ones(3), heights));
    heights(mask) = NaN;
    set(h,'ZData',heights)
    xlim([min(x_edge1(1:nbins)),max(x_edge1(1:nbins))])
    ylim([min(x_edge2(1:nbins)),max(x_edge2(1:nbins))])
    view(3)
    grid on
    xlabel([opts.fig.xlabels{1},', x_1'])
    ylabel([opts.fig.xlabels{2},', x_2'])
    zlabel('Number of Occurrences')
    
    if opts.fig.saveFiles
        fig_export(['x_hist',opts.fig.fname])
        fig_export(['x_hist',opts.fig.fname])
    end

    %% MPP 3-D histogram
    figure(numFig)
    numFig=numFig+1;
    clf
    hold on
    if isfield(opts.fig,'xlims')
        x_edge1=linspace(opts.fig.xlims.mpp(1,1),opts.fig.xlims.mpp(1,2),nbins+1);
        x_edge2=linspace(opts.fig.xlims.mpp(1,2),opts.fig.xlims.mpp(2,2),nbins+1);
    else
        x_edge1=linspace(min(min(norm_design(Umpp_final(:,1),'u',0,opts))),max(max(norm_design(Umpp_final(:,1),'u',0,opts))),nbins+1);
        x_edge2=linspace(min(min(norm_design(Umpp_final(:,2),'u',0,opts))),max(max(norm_design(Umpp_final(:,2),'u',0,opts))),nbins+1);
    end        

    x_edge1=[x_edge1(1:nbins),inf];
    x_edge2=[x_edge2(1:nbins),inf];

    edges{1}=x_edge1;
    edges{2}=x_edge2;
    hist3(norm_design(Umpp_final,'u',0,opts),'Edges',edges)
    h = get(gca,'child');
    set(h,'FaceColor','interp', 'CdataMode', 'auto');
    heights=get(h,'Zdata');
    mask=~logical(filter2(ones(3), heights));
    heights(mask) = NaN;
    set(h,'ZData',heights)

    if isfield(opts.fig,'ylims')
        zlim([opts.fig.ylims.mpp(1,1),opts.fig.ylims.mpp(1,2)])
    end
    
    if isfield(opts.fig,'caxis')
        caxis([opts.fig.caxis.mpp(1,1),opts.fig.caxis.mpp(1,2)])
    end
    
    unorm=norm_design(opts.design.udet,'u',0,opts);
    h=plot3([unorm(1),unorm(1)],[unorm(2),unorm(2)],zlim,'--ro','LineWidth',1,'MarkerFaceColor','r','MarkerSize',2);
    hh=legend(h,'$\mathbf{u}_{det}$',-1);
    set(hh,'Interpreter','latex');

    view(3)
    grid on
    xlabel('U_1^{MPP}')
    ylabel('U_2^{MPP}')
    zlabel('Number of Occurrences')
    if isfield(opts.fig,'xlims')
        xlim([opts.fig.xlims.mpp(1,1),opts.fig.xlims.mpp(1,2)])
        ylim([opts.fig.xlims.mpp(2,1),opts.fig.xlims.mpp(2,2)])
    end
    
    if opts.fig.saveFiles
        fig_export(['mpp_hist',opts.fig.fname])
        fig_export(['mpp_hist',opts.fig.fname])
        
        legend off
        grid off
        set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[],'zlabel',[])

        view([0,0,1])
        fig_export(['mpp_hist',opts.fig.fname,'_top'],0.5/3)

        view([0,-1,0])
        fig_export(['mpp_hist',opts.fig.fname,'_u1'],0.5/3)

        view([1,0,0])
        fig_export(['mpp_hist',opts.fig.fname,'_u2'],0.5/3)
    end
end

end

