Wab = -1000;
Aab=0;
L = 0.16;
b= 0.06;
Mp= 2.5;
Mbd =3;
g=9.81;

theta = [0:1:180];

Beta = asind((b.*(sind(theta)))/L) ;
Wdb = (29.27.*sind(90-theta))/(cosd(Beta)) ; 
alphaBD = ( (657.976.*cosd(90-theta)) - (0.16.*(Wdb.^2).*sind(Beta)))/(-0.16.*cosd(Beta));

Ad = ((-657.976.*(sind(90-theta))) - (0.16.*alphaBD.*sind(Beta))-(0.16.*(Wdb.^2).*cosd(Beta)));

Dy = Mp.*(Ad+g);

Dx =  -Mp.*(Ad+g).*tand(Beta);
figure(1)
plot(theta, Dy, theta, Dx);
title('Reaction Forces Dy and Dx at Theta')
xlabel('\Theta (degrees)')
ylabel('DyDx (Newtons)')
legend({'Dy', 'Dx'}, 'Location', 'southeast')
grid

Ay= Ad - (0.08.*alphaBD.*cosd(Beta)) - (0.08.*((Wdb.^2).*sin(Beta)));
Ax = (-0.08.*alphaBD).*sind(Beta) + (0.08.*(Wdb.^2).*cosd(Beta));

By = (Mbd.*Ay) + (Mbd.*g) + (Mp.*(Ad+g));
Bx = (Mbd.*Ax) - ((Dx.*((1/12).*Mbd.*(L.^2))).*Ax) + ((Mp.*(Ad+g)).*tand(Beta));

figure(2)
plot(theta, By, theta, Bx);
title('Reaction Forces By and Bx at Theta')
xlabel('\Theta (degrees)')
ylabel('ByBx (Newtons)')
legend({'By', 'Bx'}, 'Location', 'southeast')
grid
