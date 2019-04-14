function []=fig_export(fname,ratio)
% Save figure to image file
%
% Inputs:
% fname = file name
% ratio = divide figure size by ratio (e.g. 0.5 is default for side by
% side)
%
% Outputs:
% none
%==========================================================================

% Get handles to current figure and axes
fig=gcf;
ax=gca;

if nargin==1
    ratio=0.5;
end

figType='doc';
% Figure size
switch lower(figType)
    case {'doc'}
        xSize=8.255*(2*ratio);              % Figure width
        ySize=8.255*(2*ratio);              % Figure height
        fig.PaperUnits='centimeters';
        fontSize=8;
    case {'beamer'}
        xSize=5.65*(2*ratio);              % Figure width
        ySize=5.65*(2*ratio);              % Figure height
        fig.PaperUnits='centimeters';
        fontSize=6;
    otherwise
        error('Unrecognized figure type')
end

% Set renderer
opengl software

% Format figure
movegui(fig,'center')
fig.PaperPositionMode='manual';
fig.PaperPosition=[0 0 xSize ySize];
movegui(fig,'center')

% Format axes
ax.FontSize=fontSize;
ax.FontName='Arial';

% Print file
print(fig,[pwd,'\figures\',fname],'-dpng','-r300')

% Save matlab file
savefig([pwd,'\figures\',fname])

end
