function h = histfit_YQdefined(data,dist,blimit,binNum)
%HISTFIT Histogram with superimposed fitted normal density.
%   HISTFIT(DATA,NBINS) plots a histogram of the values in the vector DATA,
%   along with a normal density function with parameters estimated from the
%   data.  NBINS is the number of bars in the histogram. With one input
%   argument, NBINS is set to the square root of the number of elements in
%   DATA. 
%
ax = newplot;

data = data(:);
data(isnan(data)) = [];
n = numel(data);
pd = fitdist(data,dist);

% Find range for plotting
q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
x = linspace(q(1),q(2));
if ~pd.Support.iscontinuous
    % For discrete distribution use only integers
    x = round(x);
    x(diff(x)==0) = [];
end

% Do histogram calculations
[bincounts,binedges] = histcounts(data,binNum,'BinLimits',blimit);
bincenters = binedges(1:end-1)+diff(binedges)/2;

% Plot the histogram with no gap between bars.
hh = bar(ax, bincenters,bincounts,1);

% Normalize the density to match the total area of the histogram
binwidth = binedges(2)-binedges(1); % Finds the width of each bin
area = n * binwidth;
y = area * pdf(pd,x);

% Overlay the density
np = get(ax,'NextPlot');
set(ax,'NextPlot','add');
hh1 = plot(ax, x,y,'r-','LineWidth',2);

if nargout == 1
  h = [hh; hh1];
end

set(ax,'NextPlot',np);
ylim([0,1.5*max(bincounts,[],'all')]);