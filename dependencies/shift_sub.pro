function shift_sub, image, x0, y0, cubic=cubic
;+
; NAME: SHIFT_SUB
; PURPOSE:
;     Shift an image with subpixel accuracies
; CALLING SEQUENCE:
;      Result = shift_sub(image, x0, y0)
;-

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
 
if fix(x0)-x0 eq 0. and fix(y0)-y0 eq 0. then return, shift(image, x0, y0)

s =size(image)
x=findgen(s[1])#replicate(1., s[2])
y=replicate(1., s[1])#findgen(s[2])
x1= (x-x0)>0<(s[1]-1.)  
y1= (y-y0)>0<(s[2]-1.)

if keyword_set(cubic) then begin
	return, interpolate(image, x1, y1, cubic=cubic, /double, MISSING=median(image, /double, /even))
endif else begin
	return, interpolate(image, x1, y1, /double, MISSING=median(image, /double, /even)); bilinear interpolation by default
endelse

end