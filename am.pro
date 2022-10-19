PRO AM
x = findgen(1000*!PI, start=-500*!PI) * (1./500.)
E = 0.1 + (sin(x))^2 + 0.35*EXP(-0.5*(x/0.4)^2)
H = cos(21.5*x)
G = EXP(-0.5*((x-0.15)/0.10)^2)
y = E * H * (1-0.5*G) + (1./6.)*G
plot, x, y, THICK=7, TITLE="AM", YRANGE=[-2, 2], XRANGE=[-2.5, 2.5]
END
