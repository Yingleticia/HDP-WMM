

yy = linspace(1,2,200); % grid values
alpha = 12.0072;
beta = 63.0269;
pp = wblpdf(yy, beta^(1/alpha),alpha);
cc = wblcdf(yy, beta^(1/alpha),alpha);

figure; plot(yy,pp)
figure; plot(yy,cc)