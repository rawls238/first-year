clear all;
close all;

% Part (a)
n=5000;
a=1.4;
b=0.3;
x=zeros(n,1);
y=zeros(n,1);

for i=2:n
   x(i) = a - (x(i-1))^2 + b*y(i-1);
   y(i) = x(i-1);
end

plot(1:5000,x,'Color',[0.5,0.5,0.5])
print('Fig_a1','-dpng')

plot(x,y,'Color',[0.5,0.5,0.5])
print('Fig_a2','-dpng')

% Part (b)
histogram(x,40,'FaceColor',[0.5,0.5,0.5])
print('Fig_b','-dpng')

% Part (c)
xx=zeros(n,1);
yy=zeros(n,1);
xx(1)=10^(-8);
yy(1)=10^(-8);

for i=2:n
   xx(i) = a - (xx(i-1))^2 + b*yy(i-1);
   yy(i) = xx(i-1);
end

plot(1:5000,xx,'Color',[0.5,0.5,0.5])
print('Fig_c1','-dpng')

plot(xx,yy,'Color',[0.5,0.5,0.5])
print('Fig_c2','-dpng')

histogram(xx,40,'FaceColor',[0.5,0.5,0.5])
print('Fig_c3','-dpng')

% Part (d)
diff=log(abs(x-xx));
plot(1:5000,diff,'Color',[0.5,0.5,0.5])
print('Fig_d','-dpng')
