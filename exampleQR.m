function exampleQR(figDir)
% Show an example for a change in a distribution, and how it maps onto the QR trends
  
k = 10000;
D1 = pearsrnd(0, 2, -0.5, 3, k, 1);
D2 = pearsrnd(0, 1.5, -0.2, 3, k, 1);

P = 1:99;
cdf1 = prctile(D1, P); 
cdf2 = prctile(D2, P);

clf
plot(cdf1, P, cdf2, P, 'linewidth', 3);
set(gca,'fontsize',20)
orient landscape
print('-dpsc',[figDir '/exampleCDF.ps'])

dQ = cdf2 - cdf1;
clf
plot(P,dQ, 'k', 'linewidth',3)
set(gca,'fontsize',20)
orient landscape
print('-dpsc',[figDir '/exampleQR.ps'])

xi = linspace(-10,5, 100);
pdf1 = ksdensity(D1, xi, 'bandwidth', 0.5);
pdf2 = ksdensity(D2, xi, 'bandwidth', 0.5);
clf
plot(xi,pdf1,xi,pdf2,'linewidth',5)
set(gca,'fontsize',40)
orient landscape
print('-dpsc',[figDir '/examplePDF.ps'])