x = findgen(1000, START=-500) * !PI / 2400
k = !PI * x
;y = cos(24*k)+(0.5)*(cos(22*k)+cos(26*k))+(0.33333333)*(cos(20*k)+cos(28*k))+(0.25)* $
;	(cos(18*k) + cos(30*k))+(0.2)*(cos(16*k)+cos(32*k))+(0.16666666)*(cos(14*k)+cos(34*k))

arg = []
for i=1,11 do begin
	arg = [arg, 2 * i * !PI * x]
endfor

coeffs = [1./6., 1./5., 1./4., 1./3., 1./2., 1., 1./2., 1./3., 1./4., 1./5., 1./6.]

y=0
for i=1,11 do begin
	y += coeffs[i-1]*cos(arg[i-1])
endfor

plot, x, y

