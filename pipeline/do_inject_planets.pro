pro do_inject_planets, obj_name, cube_folder, n_planets, planet_contrast, pxscale_sx,$
	pxscale_dx, ct, do_cen_filter, planet_r=planet_r, planet_theta=planet_theta,$
	planet_x=planet_x, planet_y=planet_y, use_gauss=use_gauss, silent=silent,$
	truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx, runs=runs
	
compile_opt idl2
newline = string(10B)
	
;Do this for runs eq 1 and runs eq 2
if runs lt 3 then begin
   output_folder = cube_folder + 'processed_left'
   truenorth = truenorth_sx
   pxscale = pxscale_sx
; do this for runs 3 and 4
endif else begin
   output_folder = cube_folder + 'processed_right'
   truenorth = truenorth_dx
   pxscale = pxscale_dx
endelse
print, 'Pixel Scale:', pxscale, newline

; Do this for runs eq 1 and runs eq 3
if runs mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'

obj_cube = readfits(output_folder + dither_folder + obj_name + string(ct) +$
	'_cube_skysub_cen_clean.fits')
	
;unsaturated data can just use the pupil image
ref = readfits(output_folder + dither_folder + obj_name + string(ct) +$
	'_pupil.fits')
	
if use_gauss eq 1 then ref = gauss2dfit(ref, /tilt)

x_size = (size(ref))[1] & y_size = (size(ref))[2]
;half-sizes
x_hsize = x_size / 2.
y_hsize = y_size / 2.

restore, filename = output_folder + dither_folder + obj_name + string(ct) +$
	'_parang_clean.sav'

;derotate
for ii=0, (size(obj_cube))[3]-1 do begin
   if not keyword_set(silent) or silent eq 0 then print, 'Derotating by ',$
    	-angles[ii]-truenorth
    	
   obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii]-truenorth, /interp)
endfor

big_ref=obj_cube[*,*,0]
big_ref[*]=0

;build bigger reference image
for xx=0, (size(ref))[1] - 1 do begin
   for yy=0, (size(ref))[1] - 1 do begin

      big_ref[x_hsize-(size(ref))[1]/2.+xx,y_hsize-(size(ref))[1]/2.+yy] = ref[xx,yy]

   endfor
endfor


;inject planets
for ii=0, n_planets-1 do begin
   Print, 'Injecting planet', ii, ' in run: ', runs, newline
   big_ref_ii=big_ref*planet_contrast[ii]
      ; If angle and radius are given
      if keyword_set(planet_theta) and keyword_set(planet_r) then begin
         ;Sin and Cos take their arguments in [rad]
     	 xshift= planet_r[ii] * (1./pxscale) * Cos(planet_theta[ii])
     	 yshift= planet_r[ii] * (1./pxscale) * Sin(planet_theta[ii])
      endif else begin; If x and y are given
         xshift = planet_x-x_hsize
         yshift = planet_y-y_hsize
      endelse

      big_ref_ii=fshift(big_ref_ii, xshift,yshift)
   
   for jj=0,(size(obj_cube))[3] - 1 do obj_cube[*,*,jj] += big_ref_ii
   
endfor; ii for

;rotate back to pupil stabilized orientation
for ii=0, (size(obj_cube))[3]-1 do begin
   if not keyword_set(silent) or silent eq 0 then print, 'Rotating by ', angles[ii]
   obj_cube[*,*,ii]=rot(obj_cube[*,*,ii],angles[ii]+truenorth,/interp)
endfor

writefits, output_folder + dither_folder + obj_name + string(ct) +$
	'_cube_skysub_cen_clean_inj.fits', obj_cube

END; do_inject_planets procedure end