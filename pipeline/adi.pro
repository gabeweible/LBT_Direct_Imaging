pro adi, obj_name, half_cropped_sz, nod, cube_folder, use_injection, do_destripe, filter, ct,$
	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band
; silent will supress all the "Rotating by ..." prints.
compile_opt idl2
newline = string(10B)

adi_suffix = cube_folder + 'combined/' + obj_name + 'ct_' + string(ct)$
	+ 'filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj) + '_uncert_' +$
	string(uncert)

if not keyword_set(magnify) then magnify = 1; default to magnification
min_pxscale = min([pxscale_sx, pxscale_dx])

if band eq 'Lp' then fwhm = 9.6 ; px
if band eq 'M' then fwhm = 13.0773; px

if nod eq 'sx_only' then begin
	
	; loop through runs/nods
	for runs=1,2 do begin
	
		;call do_adi for each run/nod
		adiframe = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, do_destripe, filter, ct,$
			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
			fwhm=fwhm)
		
		; create our cube of adi nods
		if runs eq 1 then begin
			adinods = adiframe
		endif else begin
			adinods = [[[adinods]], [[adiframe]]]
		endelse
	endfor
	
	; write the nods as separate frames in one cube
	writefits, strcompress(adi_suffix + '.fits', /rem), adinods
	
	left_adi = adinods

	left_adi_mean = (left_adi[*,*,0] + left_adi[*,*,1]) / 2
	writefits, strcompress(adi_suffix + '_left_adi.fits', /rem), left_adi_mean

	; write all the sides and the total
	writefits, strcompress(adi_suffix + '_adi_nod1.fits', /r), left_adi[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod2.fits', /r), left_adi[*,*,1]


	; unsaturated data can just use a median pupil image
	ref_file = cube_folder + 'processed_left/dith1/' + obj_name + string(ct) + '_pupil.fits'

	; call find_sources, if wanted.
	if fs eq 1 then begin

		find_sources,strcompress(adi_suffix + '_left_adi.fits',/r),$
			reference=ref_file, platescale=min_pxscale,$
			correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
			fwhm=fwhm,ct=ct,filter=filter; rough PSF size in px for M???
	  
	endif; fs set if
	
endif; nod eq 'sx_only' if

; inject into all four nods and stack
if nod eq 'total' then begin
	
	; loop through runs/nods
	for runs=1,4 do begin
	
		;call do_adi for each run/nod
		adiframe = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, do_destripe, filter, ct,$
			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
			fwhm=fwhm)
		
		; create our cube of adi nods
		if runs eq 1 then begin
			adinods = adiframe
		endif else begin
			adinods = [[[adinods]], [[adiframe]]]
		endelse
	endfor
	
	; write the nods as separate frames in one cube
	writefits, strcompress(adi_suffix + '.fits', /rem), adinods
	
	; split into left and right, also
	left_adi = adinods[*,*,0:1]
	right_adi = adinods[*,*,2:3]

	; Total ADI (scale left/SX frames to match pxscale of right/DX frames:
	; Get the factor we need to magnify the side with the greater pxscale by
	; Bigger pxscale means less pixels per arcsec so we need to scale it up to get
	; more pixels!
	if magnify eq 1 then begin

		mag_factor = max([pxscale_sx, pxscale_dx]) / min([pxscale_sx, pxscale_dx])
		old_dim = (size(right_adi[*,*,0]))[1]
		new_dim = old_dim * mag_factor
		dim_diff = new_dim - old_dim
		start_i = dim_diff / 2
		end_i = (new_dim-1) - start_i

		if pxscale_sx gt pxscale_dx then begin
		
			new_xdim = (size(left_adi[*,*,0]))[1] * mag_factor
			new_ydim = (size(left_adi[*,*,0]))[2] * mag_factor
			
			new_left_adi = CONGRID(left_adi[*,*,0], new_xdim, new_ydim, /INTERP)
			
			new_left_adi = [[[left_adi]], [[CONGRID(left_adi[*,*,1],$
				new_xdim, new_ydim, /INTERP)]]]
			
			left_adi = new_left_adi[start_i:end_i, start_i:end_i, *]
				
		endif else begin
			
			new_xdim = (size(right_adi[*,*,0]))[1] * mag_factor
			new_ydim = (size(right_adi[*,*,0]))[2] * mag_factor
			
			new_right_adi = CONGRID(right_adi[*,*,0], new_xdim, new_ydim, /INTERP)
			
			new_right_adi = [[[new_right_adi]], [[CONGRID(right_adi[*,*,1],$
				new_xdim, new_ydim, /INTERP)]]]
				
			right_adi = new_right_adi[start_i:end_i, start_i:end_i, *]
				
		endelse
	
	endif; magnify if

	; get means for each side
	right_adi_mean = (right_adi[*,*,0] + right_adi[*,*,1]) / 2
	writefits, strcompress(adi_suffix +  '_right_adi.fits', /rem), right_adi_mean

	left_adi_mean = (left_adi[*,*,0] + left_adi[*,*,1]) / 2
	writefits, strcompress(adi_suffix + '_left_adi.fits', /rem), left_adi_mean

	; Average the averages for both sides
	total_adi = (left_adi_mean + right_adi_mean) / 2

	; write all the sides and the total
	writefits, strcompress(adi_suffix + '_adi_nod1.fits', /r), left_adi[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod2.fits', /r), left_adi[*,*,1]
	writefits, strcompress(adi_suffix + '_adi_nod3.fits', /r), right_adi[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod4.fits', /r), right_adi[*,*,1]

	writefits, strcompress(adi_suffix + '_total_adi.fits', /r), total_adi

	; unsaturated data can just use a median pupil image
	ref_file = cube_folder + 'processed_left/dith1/' + obj_name + string(ct) + '_pupil.fits'

	; call find_sources, if wanted.
	if fs eq 1 then begin

		find_sources,strcompress(adi_suffix + '_total_adi.fits',/r),$
			reference=ref_file, platescale=min_pxscale,$
			correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
			fwhm=fwhm,ct=ct,filter=filter
	  
	endif; fs set if
	
endif; end nod eq 'total' if

;for only one nod
if (nod eq '1') or (nod eq '2') or (nod eq '3') or (nod eq '4') then begin
	; set the run and call do_adi for our nod/run
	runs = fix(nod)
	adinods = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, do_destripe, filter, ct,$
		do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
		silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
		pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
		fwhm=fwhm)
	
	; magnify, if needed.
	if magnify eq 1 then begin
		mag_factor = max([pxscale_sx, pxscale_dx]) / min([pxscale_sx, pxscale_dx])
		old_dim = (size(adinods))[1]
		new_dim = old_dim * mag_factor
		dim_diff = new_dim - old_dim
		start_i = dim_diff / 2
		end_i = (new_dim-1) - start_i

		if (pxscale_sx gt pxscale_dx) && (fix(nod) lt 3) then begin
		
			new_xdim = (size(adinods))[1] * mag_factor
			new_ydim = (size(adinods))[2] * mag_factor
			adinods = CONGRID(adinods, new_xdim, new_ydim, /INTERP)
			adinods = adinods[start_i:end_i, start_i:end_i]
			
		endif
		
		if (pxscale_dx gt pxscale_sx) && (fix(nod) ge 3) then begin
		
			new_xdim = (size(adinods))[1] * mag_factor
			new_ydim = (size(adinods))[2] * mag_factor
			adinods = CONGRID(adinods, new_xdim, new_ydim, /INTERP)
			adinods = adinods[start_i:end_i, start_i:end_i]
			
		endif
	endif; magnify if
	
	; write the adi fits for our nod
	writefits, strcompress(adi_suffix + '_adi_nod' + nod + '.fits', /r), adinods
	
endif; end one-nod if

print, 'Finished.'
end; end of procedure
