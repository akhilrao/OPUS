function [h,X1,X2,yi] = PlotDist(xim,yim,GM,fn)
[X1,X2] = meshgrid(xim,yim);
X = [X1(:) X2(:)];

yi=zeros(length(X),1);

figure(fn)
for i = 1:length(GM)
    if ~isempty(GM(i).state)
    if isfinite(GM(i).state(1))
        yi = yi+GM(i).c*mvnpdf(X,GM(i).state,GM(i).P);
    end
    end
end
yi = reshape(yi,length(yim),length(xim));

h=surf(X1,1000*X2,yi);
xlabel('Orbit Radius [km]')
ylabel('Debris Radius [m]')
view(0,90)

set(h,'edgecolor','none')
end