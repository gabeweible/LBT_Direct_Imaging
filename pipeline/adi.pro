pro adi, obj_name, half_cropped_sz, nod, cube_folder, use_injection, filter, keep_number,$
	do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
	silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
	pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, band=band, fwhm=fwhm,$
	combine_type=combine_type, szbin=szbin, bin_type=bin_type,$
	peak_thresh=peak_thresh, stddev_thresh=stddev_thresh,$
	mask=mask, outmask=outmask, do_smooth=do_smooth
	
;;;
;;; band='Lp' or 'M' for LMIRCam, nod='total' for all data (2 nods, 2 apertures)
;;;
; silent will supress all the "Rotating by ..." prints.
compile_opt idl2
newline = string(10B)

if not keyword_set(szbin) then szbin=1
if not keyword_set(bin_type) then bin_type='median'

bin=szbin

adi_suffix = cube_folder + 'combined/' + obj_name + '_keepnumber_' + string(keep_number)$
	+ '_filt_'  +  string(filter) + '_neg_inj_' + string(neg_inj) + '_uncert_' +$
	string(uncert) + '_szbin_' + string(szbin) + '_type_'+bin_type+'_comb_type_'+$
	combine_type+'_normal_'+string(normal)+'_peak_thresh_'+string(peak_thresh)+$
	'_stddev_thresh_'+string(stddev_thresh)

if not keyword_set(magnify) then begin
	magnify = 0; default to no magnification
endif

;if band eq 'Lp' then fwhm = 9.55415 ; px
;if band eq 'M' then fwhm = 13.0773; px

if nod eq 'sx_only' then begin
	
	; loop through runs/nods
	for runs=1,2 do begin
	
		;call do_adi for each run/nod
		adiframe = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, filter, keep_number,$
			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
			fwhm=fwhm, combine_type=combine_type, szbin=szbin,$
			bin_type=bin_type, do_smooth=do_smooth)
		
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

	; left combination only (three options)
	left_adi_mean=mean(left_adi, dim=3, /double, /nan)
	writefits, strcompress(adi_suffix + '_left_adi.fits', /rem), left_adi_mean
	
	; write all the sides and the total
	writefits, strcompress(adi_suffix + '_adi_nod1.fits', /r), left_adi[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod2.fits', /r), left_adi[*,*,1]

	; unsaturated data can just use a median pupil image
	ref_file = cube_folder + 'processed_left/dith1/' + obj_name + '_' + string(keep_number) + '_pupil.fits'

	; call find_sources, if wanted.
	if fs eq 1 then begin

		if magnify eq 1 then pxscale = min_pxscale else pxscale = pxscale_sx
		
		find_sources,strcompress(adi_suffix + '_left_adi.fits',/r),$
			reference=ref_file, platescale=pxscale,$
			correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
			fwhm=fwhm,keep_number=keep_number,filter=filter; rough PSF size in px for M???
	  
	endif; fs set if
	
endif; nod eq 'sx_only' if

if nod eq 'dx_only' then begin
	
	; loop through runs/nods
	for runs=3,4 do begin
	
		;call do_adi for each run/nod
		adiframe = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, filter, keep_number,$
			do_cen_filter, coadd, adi_cube_nframes, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
			fwhm=fwhm, combine_type=combine_type, szbin=szbin,$
			bin_type=bin_type, peak_thresh=peak_thresh,$
			stddev_thresh=stddev_thresh, adi_suffix=adi_suffix, do_smooth=do_smooth)
		
		; create our cube of adi nods
		if runs eq 3 then begin
			adinods = adiframe
			numbers_of_frames = adi_cube_nframes
		endif else begin
			adinods = [[[adinods]], [[adiframe]]]
			numbers_of_frames = [[numbers_of_frames], adi_cube_nframes]
		endelse

	endfor
	
	; write the nods as separate frames in one cube
	writefits, strcompress(adi_suffix + '.fits', /rem), adinods
	
	if (mask ge 1) or (outmask le half_cropped_sz) then begin
		for xx = 0, 2.*half_cropped_sz-1 do begin
			for yy=0, 2.*half_cropped_sz-1 do begin
				test_rad = sqrt((xx - half_cropped_sz)^2. + (yy-half_cropped_sz)^2.)
				if test_rad le mask then adinods[xx,yy, *] = 0.0
				
				if test_rad ge outmask then adinods[xx,yy,*] = 0.0
			endfor
		endfor
	endif
	
	; right combination only
	right_adi_mean = (numbers_of_frames[0]*adinods[*,*,0] + numbers_of_frames[1]*adinods[*,*,1])/total(numbers_of_frames)
	
	sz = 2*half_cropped_sz
	if do_smooth gt 0 then begin
		PSF_lp = psf_Gaussian(npix=11, FWHM=[do_smooth, do_smooth],$
			/double, /normalize)
		
	  	nan_indices = where(finite(right_adi_mean) ne 1, nan_count)
		if nan_count gt 0 then right_adi_mean[nan_indices] = median(right_adi_mean, /double, /even)
		
		right_adi_mean = convolve(temporary(right_adi_mean), PSF_lp)
		
		for jj=0,(size(adinods))[3]-1 do begin
			frame_jj = adinods[*,*,jj]
			
			nan_indices = where(finite(frame_jj) ne 1, nan_count)
			if nan_count gt 0 then frame_jj[nan_indices] = median(frame_jj, /double, /even)
			
			adinods[*,*,jj] = convolve(frame_jj, PSF_lp)
		endfor
	endif

	; write all the sides and the total
	writefits, strcompress(adi_suffix + '_adi_nod3.fits', /r), adinods[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod4.fits', /r), adinods[*,*,1]
	
	writefits, strcompress(adi_suffix + '_right_adi.fits', /rem), right_adi_mean
	
	print, 'PSF Width: ', fwhm
	; 33 ~ 5 * fwhm/2. in M-band
	PSFN = psf_Gaussian(npixel=73, FWHM=[fwhm/2., fwhm/2.],$
		/normalize, /double)
	
	right_adi_mean_c = convolve(right_adi_mean, PSFN)
	writefits, strcompress(adi_suffix + '_right_adi_half_conv.fits', /rem), right_adi_mean_c

	; unsaturated data can just use a median pupil image
	if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
	
	ref_file = strcompress(cube_folder + dither_folder + obj_name + '_' + string(keep_number) + '_peak_thresh_'+string(peak_thresh)+$
			'_stddev_thresh_'+string(stddev_thresh) + '_pupil.fits', /r)
			
	; call find_sources, if wanted.
	if fs eq 1 then begin
		if magnify eq 1 then pxscale = min_pxscale else pxscale = pxscale_dx

		find_sources,strcompress(adi_suffix + '_right_adi.fits',/r),$
			reference=ref_file, platescale=pxscale,$
			correction_factor=1,$
			fwhm=fwhm,keep_number=keep_number,filter=filter; rough PSF size in px for M???
	  
	endif; fs set if
	
endif; nod eq 'sx_only' if

; inject into all four nods and stack
if nod eq 'total' then begin
	min_pxscale = min([pxscale_sx, pxscale_dx])
	
	; loop through runs/nods
	for runs=1,4 do begin
	
		;call do_adi for each run/nod
		adiframe = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, filter, keep_number,$
			do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
			silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
			pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
			fwhm=fwhm, combine_type=combine_type, szbin=szbin,$
			bin_type=bin_type, do_smooth=do_smooth)
		
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
			
			new_left_adi = [[[new_left_adi]], [[CONGRID(left_adi[*,*,1],$
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
	
	; left combination only
	left_adi_mean=mean(left_adi,dim=3,/double)
	writefits, strcompress(adi_suffix + '_left_adi.fits', /rem), left_adi_mean
	
	; right combination only
	right_adi_mean=mean(right_adi,dim=3,/double)
	writefits, strcompress(adi_suffix + '_right_adi.fits', /rem), right_adi_mean
	
	total_adi_mean=mean([[[left_adi_mean]], [[right_adi_mean]]], dim=3, /double)
	writefits, strcompress(adi_suffix + '_total_adi.fits', /r), total_adi_mean

	; write all the sides and the total
	writefits, strcompress(adi_suffix + '_adi_nod1.fits', /r), left_adi[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod2.fits', /r), left_adi[*,*,1]
	writefits, strcompress(adi_suffix + '_adi_nod3.fits', /r), right_adi[*,*,0]
	writefits, strcompress(adi_suffix + '_adi_nod4.fits', /r), right_adi[*,*,1]

	; unsaturated data can just use a median pupil image
	if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
	ref_file = cube_folder + 'processed_right/'+dither_folder + obj_name + '_' + string(keep_number) + '_pupil.fits'

	; call find_sources, if wanted.
	if fs eq 1 then begin

		find_sources,strcompress(adi_suffix + '_total_adi.fits',/r),$
			reference=ref_file, platescale=min_pxscale,$
			correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),$
			fwhm=fwhm,keep_number=keep_number,filter=filter
	  
	endif; fs set if
	
endif; end nod eq 'total' if

;for only one nod
if (nod eq '1') or (nod eq '2') or (nod eq '3') or (nod eq '4') then begin
	; set the run and call do_adi for our nod/run
	runs = fix(nod)
	adinods = do_adi(obj_name, half_cropped_sz, cube_folder, use_injection, filter, keep_number,$
		do_cen_filter, coadd, fs=fs, neg_inj=neg_inj, normal=normal, uncert=uncert,$
		silent=silent, truenorth_sx=truenorth_sx, truenorth_dx=truenorth_dx,$
		pxscale_sx=pxscale_sx, pxscale_dx=pxscale_dx, magnify=magnify, runs=runs,$
		fwhm=fwhm, combine_type=combine_type, szbin=szbin,$
		bin_type=bin_type, do_smooth=do_smooth)
	
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
