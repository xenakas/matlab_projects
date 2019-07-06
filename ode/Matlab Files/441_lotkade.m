function yprime=lotkade(t,y,flag,sn,A,B,a,b,c,d)
yprime=sn*[y(1)*(A+a*y(1)+b*y(2)); y(2)*(B+c*y(1)+d*y(2))];
