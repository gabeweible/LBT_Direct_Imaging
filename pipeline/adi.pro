pro adi, obj_name, cube_folder, use_injection, do_destripe, filter, suffix, ct,$
	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, nod
; silent will supress all the "Rotating by ..." prints.
compile_opt idl2
newline = string(10B)

if not keyword_set(magnify) then magnify = 1; default to magnification
min_pxscale = min([pxscale_sx, pxscale_dx])

; inject into all four nods and stack
if nod eq 'total' then begin
	
	; loop through runs/nods
	for runs=1,4 do begin
	
				; Do this for runs eq 1 and runs eq 3
		if runs mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
		;Do this for runs eq 1 and runs eq 2
		if runs lt 3 then begin
			output_folder = cube_folder + 'processed_left'
			truenorth = truenorth_sx
		endif else begin
			truenorth = truenorth_dx
			output_folder = cube_folder + 'processed_right'
		endelse
		;call do_adi for each run/nod
		adiframe = do_adi(obj_name, cube_folder, use_injection, do_destripe, filter, suffix, ct,$
			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, nod=nod, runs=runs)
		
		; create our cube of adi nods
		if runs eq 1 then begin
			adinods = adiframe
		endif else begin
			adinods = [[[adinods]], [[adiframe]]]
		endelse
	endfor
	
	; write the nods as separate frames in one cube
	writefits, strcompress(suffix + '.fits', /rem), adinods
	
	; split into left and right, also
	left_adi = adinods[*,*,0:1]
	right_adi = adinods[*,*,2:3]

	; Total ADI (scale left/SX frames to match pxscale of right/DX frames:
	; Get the factor we need to magnify the side with the greater pxscale by
	; Bigger pxscale means less pixels per arcsec so we need to scale it up to get
	; more pixels!
	if magnify eq 1 then begin

		mag_factor = max([pxscale_sx, pxscale_dx]) / min([pxscale_sx, pxscale_dx])

		if pxscale_sx gt pxscale_dx then begin
			for i=0,1 do begin
				left_adi[*,*,i] = rot(left_adi[*,*,i], 0, mag_factor, CUBIC=-0.5)
			endfor
		endif else begin
			for j=0,1 do begin
				right_adi[*,*,j] = rot(right_adi[*,*,j], 0, mag_factor, CUBIC=-0.5)
			endfor
		endelse
	
	endif; magnify if

	; get means for each side
	right_adi_mean = mean(right_adi, dim=3)
	writefits, strcompress(suffix +  '_right_adi.fits', /rem), right_adi_mean

	left_adi_mean = mean(left_adi, dim=3)
	writefits, strcompress(suffix + '_left_adi.fits', /rem), left_adi_mean

	; Average the averages for both sides
	total_adi = (left_adi_mean + right_adi_mean) / 2

	; write all the sides and the total
	writefits, strcompress(suffix + '_adi_nod1.fits', /r), left_adi[*,*,0]
	writefits, strcompress(suffix + '_adi_nod2.fits', /r), left_adi[*,*,1]
	writefits, strcompress(suffix + '_adi_nod3.fits', /r), right_adi[*,*,0]
	writefits, strcompress(suffix + '_adi_nod4.fits', /r), right_adi[*,*,1]

	writefits, strcompress(suffix + '_total_adi.fits', /r), total_adi

	; unsaturated data can just use the pupil image
	ref_file = output_folder + dither_folder + obj_name + string(ct) + '_pupil.fits'

	; call find_sources, if wanted.
	if fs eq 1 then begin

		find_sources,strcompress(suffix + '_total_adi.fits',/r),$
			reference=ref_file, platescale=min_pxscale,$
			correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
			fwhm=8.7,ct=ct,filter=filter
	  
	endif; fs set if
	
endif else begin; else for only one nod
	; set the run and call do_adi for our nod/run
	runs = fix(nod)
	adinods = do_adi(obj_name, cube_folder, use_injection, do_destripe, filter, suffix, ct,$
		do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
		silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
		pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, nod=nod, runs=runs)
	
	; magnify, if needed.
	if magnify eq 1 then begin
		mag_factor = max([pxscale_sx, pxscale_dx]) / min([pxscale_sx, pxscale_dx])

		if (pxscale_sx gt pxscale_dx) && (fix(nod) lt 3) then begin
			adinods[*,*,0] = rot(adinods[*,*,0], 0, mag_factor, CUBIC=-0.5)
		endif
		
		if (pxscale_dx gt pxscale_sx) && (fix(nod) ge 3) then begin
			adinods[*,*,0] = rot(adinods[*,*,0], 0, mag_factor, CUBIC=-0.5)
		endif
	endif; magnify if
	
	; write the adi fits for our nod
	writefits, strcompress(suffix + '_adi_nod' + nod + '.fits', /r), adinods
	
endelse; end one nod else

print, 'Finished.'
end; end of procedure
