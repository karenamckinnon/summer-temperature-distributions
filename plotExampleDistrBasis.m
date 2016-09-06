% Plot the distributions used for the basis functions
% Fig S4

figDir = '';
  
k = 1e6;
 
% first basis: shift normal by 1
D00 = normrnd(0, 1, k, 1);
D01 = normrnd(1, 1, k, 1);
 
D10 = normrnd(0,  1, k, 1);
D11 = normrnd(0, 1.5, k, 1);

D20 = pearsrnd(0, 1, -0.5, 3, k, 1); % change in skewness of 1
D21 = pearsrnd(0, 1, 0.5, 3, k, 1);
 
D30 = pearsrnd(0, 1, 0, 2, k, 1); % change in kurtosis of 1
D31 = pearsrnd(0, 1, 0, 4, k, 1);

clf
xi = linspace(-5, 5, 101);
for ct = 1:4
	subplot(2,2,ct);
	hold on
	switch ct
		case 1
			D0 = D00;
			D1 = D01;
			ylabel('Probability density')
			text(-4.5, 0.4, '(a)', 'fontsize', 10)
		case 2
			D0 = D10;
			D1 = D11;
			text(-4.5, 0.4, '(b)', 'fontsize', 10)
		case 3
			D0 = D20;
			D1 = D21;
			ylabel('Probability density')
			xlabel('Value')
			text(-4.5, 0.4, '(c)', 'fontsize', 10)
		otherwise
			D0 = D30;
			D1 = D31;
			xlabel('Value')
			text(-4.5, 0.4, '(d)', 'fontsize', 10)
	end
	
	y0 = ksdensity(D0, xi, 'bandwidth', 0.05);
	y1 = ksdensity(D1, xi, 'bandwidth', 0.05);
	ylim([0 0.45])
	plot(xi, y0, 'color', 0.7*[1 1 1], 'linewidth', 2);
	hold on
	plot(xi, y1, 'color', 'k', 'linewidth', 2);
	set(gca, 'box', 'on')
	set(gca, 'layer', 'top')
	set(gca,'fontsize',12)
	set(gca,'fontsize',12)
	

end

orient landscape
set(gcf,'color','w')

export_fig([figDir '/exampleDistr.png'],'-m3','-a1')
