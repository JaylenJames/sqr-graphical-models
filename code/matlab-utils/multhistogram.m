function [h,ax,BigAx,hhist,pax] = multhistogram(opts,varargin)
%PLOTMATRIX Scatter plot matrix.
%   PLOTMATRIX(X,Y) scatter plots the columns of X against the columns
%   of Y.  If X is P-by-M and Y is P-by-N, PLOTMATRIX will produce a
%   N-by-M matrix of axes. PLOTMATRIX(Y) is the same as PLOTMATRIX(Y,Y)
%   except that the diagonal will be replaced by HISTOGRAM(Y(:,i)).
%
%   PLOTMATRIX(...,'LineSpec') uses the given line specification in the
%   string 'LineSpec'; '.' is the default (see PLOT for possibilities).
%
%   PLOTMATRIX(AX,...) uses AX as the BigAx instead of GCA.
%
%   [H,AX,BigAx,P,PAx] = PLOTMATRIX(...) returns a matrix of handles
%   to the objects created in H, a matrix of handles to the individual
%   subaxes in AX, a handle to big (invisible) axes that frame the
%   subaxes in BigAx, a vector of handles for the histogram plots in
%   P, and a vector of handles for invisible axes that control the
%   histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
%   that a subsequent TITLE, XLABEL, or YLABEL will be centered with
%   respect to the matrix of axes.
%
%   Example:
%       x = randn(50,3); y = x*[-1 2 1;2 0 1;1 -2 3;]';
%       plotmatrix(y)

%   Copyright 1984-2015 The MathWorks, Inc.

% Parse possible Axes input
assert(isstruct(opts), 'Opts should be structure');

[cax,args,nargs] = axescheck(varargin{:});
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 3
    error(message('MATLAB:narginchk:tooManyInputs'));
end
nin = nargs;

sym = '.'; % Default scatter plot symbol.
dohist = 0;

if ischar(args{nin}),
    sym = args{nin};
    [~,~,~,msg] = colstyle(sym);
    if ~isempty(msg), error(msg); end
    nin = nin - 1;
end

if nin==1, % plotmatrix(y)
    x = args{1}; y = args{1};
    dohist = 1;
elseif nin==2, % plotmatrix(x,y)
    x = args{1}; y = args{2};
else
    error(message('MATLAB:plotmatrix:InvalidLineSpec'));
end

% Post processing
if(~isfield(opts,'labels'))
    opts.labels = cellfun(@num2str, num2cell(1:size(x,2)), 'UniformOutput', false);
end
if(~isfield(opts,'xlabel'))
    opts.xlabel = 'Pairwise Maximum Mean Discrepancy';
end
x = [NaN(size(x,1),1), x];
opts.labels = [{''},opts.labels];
cols = 1;
rows = size(x,2);

%h = newplot();
%set(h,'Visible','off','color','none');

% Don't plot anything if either x or y is empty
hhist = gobjects(0);
pax = gobjects(0);
if isempty(rows) || isempty(cols),
    if nargout>0, h = gobjects(0); ax = gobjects(0); BigAx = gobjects(0); end
    return
end

if ~ismatrix(x) || ~ismatrix(y)
    error(message('MATLAB:plotmatrix:InvalidXYMatrices'))
end
if size(x,1)~=size(y,1) || size(x,3)~=size(y,3),
    error(message('MATLAB:plotmatrix:XYSizeMismatch'));
end

% Create/find BigAx and make it invisible
BigAx = newplot(cax);
fig = ancestor(BigAx,'figure');
hold_state = ishold(BigAx);
set(BigAx,'Visible','off','color','none')

if any(sym=='.'),
    units = get(BigAx,'units');
    set(BigAx,'units','pixels');
    pos = get(BigAx,'Position');
    set(BigAx,'units',units);
    markersize = max(1,min(15,round(15*min(pos(3:4))/max(1,size(x,1))/max(rows,cols))));
else
    markersize = get(0,'DefaultLineMarkerSize');
end

% Create and plot into axes
ax = gobjects(rows,cols);
pos = get(BigAx,'Position');
inset = 0.17;
pos(1) = pos(1)+pos(3)*inset;
pos(3) = pos(3)*(1-inset);
pos(3) = pos(3)*1.1;
set(BigAx,'Position',pos);

width = pos(3)/cols;
height = pos(4)/rows;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];
m = size(y,1);
k = size(y,3);
xlim = zeros([rows cols 2]);
ylim = zeros([rows cols 2]);
BigAxHV = get(BigAx,'HandleVisibility');
BigAxParent = get(BigAx,'Parent');
paxes = findobj(fig,'Type','axes','tag','PlotMatrixScatterAx');


%% Get global xlim and ylim
if(isfield(opts,'maxX'))
    maxXlim = opts.maxX;
else
    maxXlim = quantile(x(:),1);
end

histogramVarargin = {max(floor(size(x,1)/5),3),'Normalization','pdf','EdgeColor','none','FaceColor',0.1*ones(1,3)};
if(size(x,1) > 10^2)
    histogramVarargin(1) = []; % Auto choose bin size
end
%{
maxYlim = 0;
for i=rows:-1:1,
    for j=cols:-1:1,
        if(j == 1)
            xCur = reshape(x(:,i,:),[m k]);
            H = histogram(xCur,histogramVarargin{:});
            maxYlim = max([maxYlim, get(H,'Values')]);
        end
    end
end
cla;
set(gca, 'visible','off');
%}

for i=rows:-1:1,
    for j=cols:-1:1,
        axPos = [pos(1) pos(2)+(rows-i)*height ...
            pos(3) height*(1-space)];
        findax = findaxpos(paxes, axPos);
        if isempty(findax),
            ax(i,j) = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',BigAxParent);
            set(ax(i,j),'visible','on');
        else
            ax(i,j) = findax(1);
        end
        %hh(i,j,:) = plot(reshape(x(:,j,:),[m k]), ...
        %    reshape(y(:,i,:),[m k]),sym,'parent',ax(i,j))';
        if(j == 1); 
            xCur = reshape(x(:,i,:),[m k]);
            hh(i,j,:) = histogram(xCur,...
                histogramVarargin{:},...
                'parent',ax(i,j)...
                )';
            if any(~isnan(xCur))
                temp = get(hh(i,j,:),'Values');
                maxYlim = max(temp(2:end))*1.1;
            else
                maxYlim = 1;
            end
            hold on;
            LW = 2;
            plot(nanmean(xCur)*ones(2,1), [0,maxYlim],'-r',...
                'LineWidth', LW,...
                'parent', ax(i,j),...
                'color','r');
            plot(quantile(xCur,0)*ones(2,1), [0,maxYlim],':b',...
                'LineWidth', LW,...
                'parent', ax(i,j),...
                'color','b');
            plot(quantile(xCur,0.5)*ones(2,1), [0,maxYlim],'--b',...
                'LineWidth', LW,...
                'parent', ax(i,j),...
                'color','b');
            plot(quantile(xCur,1)*ones(2,1), [0,maxYlim],':b',...
                'LineWidth', LW,...
                'parent', ax(i,j),...
                'color','b');
            hold off;
            
            % Set global axes limits
            set(ax(i,j),'xlim',[0,maxXlim]);
            if(~isnan(maxYlim))
                set(ax(i,j),'ylim',[0,maxYlim]);
            end
            
            % Set ylabel and xlabel
            labelFontSize = 12;
            set(ax(i,j), 'FontSize',labelFontSize-2);
            
            yLab = ylabel(ax(i,j), opts.labels{i});
            ylp = get(yLab,'Position');
            ylp(1) = -0.005;
            set(yLab,'FontSize', labelFontSize, 'Position', ylp, 'Rotation',0,'HorizontalAlignment', 'right', 'VerticalAlignment','middle');
            if(i == rows)
                xlabel(opts.xlabel,'FontSize',labelFontSize);
            end
            if(i == 1)
                legend({'Histogram', 'Mean','Min', 'Median', 'Max'},...
                    'FontSize',labelFontSize-2,'Location','SouthOutside','Orientation','Horizontal');
                %pause(0.001); % Needed to render legend
                pause(1);
                cla(ax(i,j));
                set(ax(i,j),'Visible','off');
            end
        end

    end
end


set(ax(1:rows-1,:),'xticklabel','')
set(ax(:),'yticklabel','','ytick',[])
%set(ax(:),'yminorgrid', 'on');

set(BigAx,'XTick',get(ax(rows,1),'xtick'),'YTick',get(ax(rows,1),'ytick'), ...
    'userdata',ax,'tag','PlotMatrixBigAx')
set(ax,'tag','PlotMatrixScatterAx');


% Make BigAx the CurrentAxes
set(fig,'CurrentAx',BigAx)
if ~hold_state,
    set(fig,'NextPlot','replacechildren')
end

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
    'String','','Visible','on')

if nargout~=0,
    h = hh;
end



function findax = findaxpos(ax, axpos)
tol = eps;
findax = [];
for i = 1:length(ax)
    axipos = get(ax(i),'Position');
    diffpos = axipos - axpos;
    if (max(max(abs(diffpos))) < tol)
        findax = ax(i);
        break;
    end
end
