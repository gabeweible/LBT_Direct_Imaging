pro inject_planets, obj_name, cube_folder, n_planets, planet_contrast, pxscale_sx,$
	pxscale_dx, ct, do_cen_filter, planet_r=planet_r, planet_theta=planet_theta,$
	planet_x=planet_x, planet_y=planet_y, use_gauss=use_gauss, silent=silent,$
	truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, nod=nod
	
;'HII1348', '~/OneDrive/Research/HII1348/testing'
compile_opt idl2
newline = string(10B)

print, 'Number of Planets:', n_planets

if keyword_set(planet_r) and keyword_set(planet_theta) then begin
	print, 'Planet Radii:', planet_r
	print, 'Planet Thetas:', planet_theta
endif else begin
	if keyword_set(planet_x) and keyword_set(planet_y) then begin
		print, 'Planet x:', planet_x
		print, 'Planet y:', planet_y
	endif else begin
		print, 'Error: no radii and thetas OR x and y values for the injection.' & stop
	endelse
endelse

print, 'Planet Contrast:', planet_contrast

if nod eq 'total' then begin
	for runs=1,4 do begin
	
		do_inject_planets, obj_name, cube_folder, n_planets, planet_contrast, pxscale_sx,$
			pxscale_dx, ct, do_cen_filter, planet_r=planet_r, planet_theta=planet_theta,$
			planet_x=planet_x, planet_y=planet_y, use_gauss=use_gauss, silent=silent,$
			truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, runs=runs
			
	endfor; runs for
endif else begin

	runs = fix(nod)
	do_inject_planets, obj_name, cube_folder, n_planets, planet_contrast, pxscale_sx,$
			pxscale_dx, ct, do_cen_filter, planet_r=planet_r, planet_theta=planet_theta,$
			planet_x=planet_x, planet_y=planet_y, use_gauss=use_gauss, silent=silent,$
			truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, runs=runs
			
endelse

end; inject_planets procedure end
