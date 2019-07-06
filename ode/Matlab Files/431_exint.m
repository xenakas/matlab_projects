%% Solution of Exercise 4.1(iii)

f=inline('exp(-(t-s)).*sin(s)','t','s');

for j=0:100;
     t(j+1)=7*j/100;
     T(j+1)=quad(f,0,t(j+1),[],[],t(j+1));
end

plot(t,T)
