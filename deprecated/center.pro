pro center, obj_name, stripe, cube_folder, do_block_right, do_block_left,$
	do_block_bottom, do_block_top, do_cen_filter, do_cen_median_sub, half_cropped_sz, aperture,$
	ct,good_framei=good_framei, align_radius_thresh=align_radius_thresh
; 'HII1348', '~/OneDrive/Research/HII1348/testing', 1, 1, 1, 1, 1, 0 for current run (let's try cen_median_sub later)
compile_opt idl2

; TEMPORARY!!! ADD KWARGS!!!
	align_filter=0 ;high-pass filter 
		align_filter_width=20.
	align_smooth=0
		align_smooth_width=5.
		align_median_sub=do_cen_median_sub
		do_destripe=0

if not keyword_set(align_radius_thresh) then align_radius_thresh=50

; I think this is the second time that the runs matter (we are now post-split)
; HII 1348 CENTER 1024 PIXEL ROWS (DOUBLE-SIDED)
if (stripe eq 'center1024') and (aperture eq 'both') then for runs=1,4 do begin
   ; Do this for runs eq 1 and runs eq 3
   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
   
   ;Do this for runs eq 1 and runs eq 2
   if runs lt 3 then output_folder = cube_folder + '/processed_left/' else output_folder = cube_folder + '/processed_right/'

   print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
   print, 'Reading in a sky_sub cube post-split...'
   obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
   print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
   restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
   print, 'Post-split flags and angles read! Changing some stuff in the cube...'
   obj_cube[where(finite(obj_cube) eq 0)] = 0.
   print, 'Some stuff in the cube has been changed! Creating the obj_cube copy...'
   
   print, 'Duplicating cube...'
	obj_cube2=obj_cube
	
	if do_destripe eq 1 then begin
	print, 'destriping 90 degrees...'
	for ii=0, (size(obj_cube2))[3]-1 do obj_cube2[*,*,ii]=destripe(obj_cube2[*,*,ii],90.,clip_level=0.0,/nodisp)
	print, 'destriping 0 degrees...'
	for ii=0, (size(obj_cube2))[3]-1 do obj_cube2[*,*,ii]=destripe(obj_cube2[*,*,ii],0.,clip_level=0.0,/nodisp)
	endif
	
	obj_cube2=obj_cube2 > 0. ;truncate negative values to not get confused on the negative PSF
	
	;block edges
	;obj_cube2[0:half_cropped_sz-50,*,*]=obj_cube2[0:20,*,*]/1000000.
	;obj_cube2[*,*,,*]=obj_cube2[*,0:20,*]/1000000.
	;obj_cube2[*,107:127,*]=obj_cube2[*,107:127,*]/1000000.
	;obj_cube2[107:127,*,*]=obj_cube2[107:127,*,*]/1000000.
	
	;print,'Reducing negatives...'
	;if 1 then begin
	;for ii=0,(size(obj_cube2))[3]-1 do begin
	;	frame=obj_cube2[*,*,ii]
	;		frame[where(frame lt 0)]=frame[where(frame lt 0)]/100000000000000. ;reduce negative values so CC doesn't get confused on subtracted PSF.
	;	obj_cube2[*,*,ii]=frame
	;endfor
	;endif
	
	
	w=half_cropped_sz
	
	obj_cube2=obj_cube2[half_cropped_sz-w:half_cropped_sz+w-1,half_cropped_sz-w:half_cropped_sz+w-1,*]
	
	half_sizex=w;half_cropped_sz
	half_sizey=w;half_cropped_sz
	
	if align_filter eq 1 then begin
	Print,'Filtering frames...'
	for ii=0, (size(obj_cube))[3]-1 do begin
		frame=obj_cube2[*,*,ii]
		frame=frame-smooth(frame,align_filter_width,/NAN)
		obj_cube2[*,*,ii]=frame
	endfor
	endif
	if align_smooth eq 1 then  begin
	Print,'Smoothing frames...'
	for ii=0, (size(obj_cube))[3]-1 do begin
		frame=obj_cube2[*,*,ii]
		frame=smooth(frame,align_smooth_width,/NAN)
		obj_cube2[*,*,ii]=frame
	endfor
	endif
	;if align_convolve eq 1 then  begin
	;	Print,'Convolving frames...'
	;	sz=2.0*w
	;	width=(10.*1E-6) / (8.4) * 206265. / pxscale
	;	;print, 'PSF Width: ',width
	;	PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	;	PSFN = PSF/MAX(PSF,/nan)
	;	for ii=0, (size(obj_cube))(3)-1 do begin
	;		frame=obj_cube2[*,*,ii]
	;		frame= convolve(frame, PSFN)
	;		obj_cube2[*,*,ii]=frame
	;	endfor
	;endif
		
	if align_median_sub eq 1 then medarr,obj_cube2,medframe
	for ii=0, (size(obj_cube))[3]-1 do if align_median_sub eq 1 then obj_cube2[*,*,ii]=obj_cube2[*,*,ii]-medframe
	
	
	writefits,'~/Desktop/align_testmed_unblocked.fits',median(obj_cube2,dim=3,/even)
	
	;focus only on inner part
	for xx=0,2.0*w-1 do for yy=0,2.0*w-1 do if sqrt((xx-w)^2.0+(yy-w)^2.0) gt align_radius_thresh then obj_cube2[xx,yy,*]=0.
	
	medarr,obj_cube2,cenmedian
	
	cenmedian=cenmedian>0 ;truncate negatives
	
	writefits,'~/Desktop/align_testmed.fits',cenmedian
	
	cenmedian=gauss2dfit(cenmedian)
	m=max(cenmedian,i)
	ix=i mod (2.0*w)
	iy=i/(2.0*w)
	cenmedian=shift_sub(cenmedian,w-ix,w-iy)
	
	writefits,'~/Desktop/align_testmed_gauss.fits',cenmedian
	
	cenmedian = readfits(output_folder + dither_folder + obj_name + string(ct) + '_pupil.fits')
	
	for ii=0, (size(obj_cube))[3]-1 do begin
		
		corr=crosscorr(cenmedian,obj_cube2[*,*,ii],pmax, /cut)
		if pmax[0] ge 0. then dx = half_sizex-pmax[0] else dx = -half_sizex+abs(pmax[0])
		if pmax[1] ge 0. then dy = half_sizey-pmax[1] else dy = -half_sizey+abs(pmax[1])
	
		print, 'Shifting frame', ii, ' / ', (size(obj_cube))[3]-1, ' by: ',String(-dx),String(-dy)
		obj_cube[*,*,ii]=shift_sub(obj_cube[*,*,ii],-dx,-dy)
	endfor
	
	;first find the peak of the median and shift to there
	medarr,obj_cube,medframe
	
	;if jjj eq 1 then writefits,strcompress(reduced_dir+name+'_'+side+'_cube1_align_median_PSF.fits',/rem), medframe
	medframe=medframe-smooth(medframe,25,/nan)
	for xx=0.,2.0*half_cropped_sz-1 do begin
		for yy=0.,2.0*half_cropped_sz-1 do begin
			xd=abs(xx-half_cropped_sz)
			yd=abs(yy-half_cropped_sz)
			if sqrt(xd*xd + yd*yd) gt 50 then medframe[xx,yy]=0.
		endfor
	endfor
	writefits,'~/Desktop/centest.fits',medframe
	m=max(medframe,i,/nan)
	ix=i mod (2.0*half_cropped_sz)
	iy=i/(2.0*half_cropped_sz)
	print, 'Center guess = ',ix,iy
	cntrd,medframe,ix,iy,ix,iy,7
	print,'Found center at: ',ix,iy
	;hak
	dx=float(half_cropped_sz)-ix-0.5
	dy=float(half_cropped_sz)-iy-0.5
	for ii=0,(size(obj_cube))[3]-1 do obj_cube[*,*,ii]=shift_sub(obj_cube[*,*,ii],dx,dy)


   ;we don't need all of the frames, so select some at random
	;random_inds=round((n_elements(angles))*randomu(Seed,100))

	;obj_cube2=obj_cube[*,*,random_inds]
	;angles2=angles[random_inds]


	;;first find the peak of the median and shift to there
	;medarr,obj_cube2,medframe

	;;if jjj eq 1 then writefits,strcompress(reduced_dir+name+'_'+side+'_cube1_align_median_PSF.fits',/rem), medframe
	;medframe=medframe-smooth(medframe,25,/nan)
	;for xx=0.,2.0*half_cropped_sz-1 do begin
	;	for yy=0.,2.0*half_cropped_sz-1 do begin
	;		xd=abs(xx-half_cropped_sz)
	;		yd=abs(yy-half_cropped_sz)
	;		if sqrt(xd*xd + yd*yd) gt 50 then medframe[xx,yy]=0.
	;	endfor
	;endfor
	
	;writefits,'~/Desktop/centest.fits',medframe
	;m=max(medframe,i,/nan)
	;ix=i mod (2.0*half_cropped_sz)
	;iy=i/(2.0*half_cropped_sz)
	;print, 'Center guess = ',ix,iy
	;cntrd,medframe,ix,iy,ix,iy,7
	;print,'Found center at: ',ix,iy

	;dx=float(half_cropped_sz)-ix-0.5
	;dy=float(half_cropped_sz)-iy-0.5
	;for ii=0,(size(obj_cube))[3]-1 do obj_cube[*,*,ii]=shift_sub(obj_cube[*,*,ii],dx,dy)
   ;correcting bad frames
   ;if output_folder eq '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed-left/' then begin
   ;if dither_folder eq '/dith1/' and night eq 2 then obj_cube[*, *, 36] = !values.f_nan
   ;if dither_folder eq '/dith2/' and night eq 2 then obj_cube[*, *, 24:38] = !values.f_nan
   ;endif
   
   ;bad_pixels = [37, 90, 91, 92, 93, 94] - 1
   ;goods = indgen(n_frames)
   
   ;remove, bad_pixels, goods
   ;print, angles
   ;obj_cube = obj_cube[*, *, goods]
   ;angles = angles[goods]
   ;print, angles
   
   ;endif
   
   ;if do_median_subtraction then for ii=0, n_frames-1 do obj_cube[*,*,ii] = obj_cube[*,*,ii] - medframe
   
   ;initialize a new cube to modify with the same data
   ;obj_cube2 = obj_cube
   ;print, 'Obj_cube copy created! Blocking out what needs blocked...'
   
   ;; Divide into oblivion if we're blocking a side out
   ;if do_block_right then obj_cube2[fix(0.65*2*half_cropped_sz):2*half_cropped_sz-1, *, *] *= 0.000001
   ;if do_block_left then obj_cube2[0:fix(0.35*2*half_cropped_sz), *, *] *= 0.000001
   ;;if do_block_left then obj_cube2[700:1023, *, *] *= 0.000001
   ;;if do_block_right then obj_cube2[0:450, *, *] *= 0.000001
   ;if do_block_bottom then obj_cube2[*, 0:fix(0.35*2*half_cropped_sz), *] *= 0.000001
   ;if do_block_top then obj_cube2[*,fix(0.65*2*half_cropped_sz):2*half_cropped_sz-1, *] *= 0.000001
   ;print, 'Requested blocking complete! Changing some stuff in the second cube...'
   
   ;obj_cube2[where(finite(obj_cube) eq 0)] = 0.

   ;obj_cube2[where(obj_cube2 lt 0)] *= 0.000000001 ;reduce negative values so CC doesn't get confused on subtracted PSF.
   ;   print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
   ;half_sizex = half_cropped_sz
   ;half_sizey = half_cropped_sz
   
   ;n_frames = (size(obj_cube))[3]
   ;if do_cen_filter then begin
    ;  for ii=0, n_frames - 1 do begin
     ;    print, 'Cen_filter on frame:' + string(ii)
      ;   print, 'Run:' + string(runs)
       ;  obj_cube2[*,*,ii] -= smooth(obj_cube2[*,*,ii], 17., /NAN)
        ; obj_cube2[*,*,ii] = smooth(obj_cube2[*,*,ii], 5., /NAN)
      ;endfor
   ;endif; do_cen_filter if
   
   ;print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
   ;; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
   ;;Commented out for now, let's try without this first.
   ;if do_cen_median_sub then medarr, obj_cube2, medframe
   
   ;for xx=0,2.0*half_sizex-1 do for yy=0,2.0*half_sizey-1 do if sqrt((xx-half_sizex)^2.0+(yy-half_sizey)^2.0) gt align_radius_thresh then obj_cube2[xx,yy,*]=0.
   
   	;; take a median and get rid of negatives
	;medarr, obj_cube2, cenmedian
	;cenmedian=cenmedian>0 ;truncate negatives

	;cenmedian=gauss2dfit(cenmedian); fit a gaussian to the median
	;m=max(cenmedian,i); 
	;ix=i mod (2.0*half_sizex)
	;iy=i/(2.0*half_sizey)
	;cenmedian=shift_sub(cenmedian,half_sizex-ix,half_sizey-iy); shift to center
   
   ;for ii=0, n_frames-1 do begin
      
    ;  ;Commented out for now, let's try without this first.
      ;if do_cen_median_sub then obj_cube2[*,*,ii] -= medframe
     ; 
      ;; must supply reference frame
      ;corr = crosscorr(cenmedian, obj_cube2[*,*,ii], pmax, /cut)
      ;if pmax[0] ge 0. then dx = half_sizex - pmax[0] else dx = -half_sizex + abs(pmax[0])
      ;if pmax[1] ge 0. then dy = half_sizey - pmax[1] else dy = -half_sizey + abs(pmax[1])
  ; 
      ;print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
      ;print, 'Run:' + string(runs)
      ;obj_cube[*,*,ii] = shift_sub(obj_cube[*,*,ii], -dx, -dy); shift_sub shifts frames by non-int values, which we're using to center
   ;endfor; ii for
   
   ;Write our centered FITS file
   writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', obj_cube
   print, 'Centered FITS written for run:' + string(runs)
  
endfor; runs = 1,4 for (HII 1348)


;ALCOR SECOND SET OF 512 PIXEL ROWS (SINGLE-SIDED SX OBSERVATIONS)
if stripe eq 'second512' then begin
	
	if aperture eq 'left' then for runs=1,2 do begin
		; Do this for runs eq 1 and runs eq 3
		if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
		
		output_folder = cube_folder + '/processed_left/'
		
		print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
		print, 'Reading in a sky_sub cube post-split...'
		obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
		print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
		restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
		print, 'Post-split flags and angles read! Changing some stuff in the cube...'
		obj_cube[where(finite(obj_cube) eq 0)] = 0.
		print, 'Some stuff in the cube has been changed! Creating the obj_cube copy...'
		
		;initialize a new cube to modify with the same data
		obj_cube2 = obj_cube
		print, 'Obj_cube copy created! Blocking out what needs blocked...'
		
		; Divide into oblivion if we're blocking a side out
		if do_block_right then obj_cube2[fix(0.65*2*half_cropped_sz):2*half_cropped_sz-1, *, *] *= 0.000000000001
		if do_block_left then obj_cube2[0:fix(0.35*2*half_cropped_sz), *, *] *= 0.0000000000001
		;if do_block_left then obj_cube2[700:1023, *, *] *= 0.000001
		;if do_block_right then obj_cube2[0:450, *, *] *= 0.000001
		if do_block_bottom then obj_cube2[*, 0:fix(0.35*2*half_cropped_sz), *] *= 0.000000000001
		if do_block_top then obj_cube2[*,fix(0.65*2*half_cropped_sz):2*half_cropped_sz-1, *] *= 0.000000000001
		print, 'Requested blocking complete! Changing some stuff in the second cube...'
		
		obj_cube2[where(finite(obj_cube) eq 0)] = 0.
		
		obj_cube2[where(obj_cube2 lt 0)] *= 0.000000000000001 ;reduce negative values so CC doesn't get confused on subtracted PSF.
		  print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
		half_sizex = half_cropped_sz
		half_sizey = half_cropped_sz
		
		n_frames = (size(obj_cube))[3]
		if do_cen_filter then begin
		  for ii=0, n_frames - 1 do begin
			 print, 'Cen_filter on frame:' + string(ii)
			 print, 'Run:' + string(runs)
			 obj_cube2[*,*,ii] -= smooth(obj_cube2[*,*,ii], 17., /NAN)
			 obj_cube2[*,*,ii] = smooth(obj_cube2[*,*,ii], 5., /NAN)
		  endfor
		endif; do_cen_filter if
		
		print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
		; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
		;Commented out for now, let's try without this first.
		if do_cen_median_sub then medarr, obj_cube2, medframe
		
		for xx=0,2.0*half_sizex-1 do for yy=0,2.0*half_sizey-1 do if sqrt((xx-half_sizex)^2.0+(yy-half_sizey)^2.0) gt align_radius_thresh then obj_cube2[xx,yy,*]=0.
   
		; take a median and get rid of negatives
		medarr, obj_cube2, cenmedian
		cenmedian=cenmedian>0 ;truncate negatives

		cenmedian=gauss2dfit(cenmedian); fit a gaussian to the median
		m=max(cenmedian,i); 
		ix=i mod (2.0*half_sizex)
		iy=i/(2.0*half_sizey)
		cenmedian=shift_sub(cenmedian,half_sizex-ix,half_sizey-iy); shift to center
		
		for ii=0, n_frames-1 do begin
		  
		  ;Commented out for now, let's try without this first.
		  if do_cen_median_sub then obj_cube2[*,*,ii] -= medframe
		  
		  corr = crosscorr(obj_cube2[*,*,cenmedian], obj_cube2[*,*,ii], pmax, /cut)
		  if pmax[0] ge 0. then dx = half_sizex - pmax[0] else dx = -half_sizex + abs(pmax[0])
		  if pmax[1] ge 0. then dy = half_sizey - pmax[1] else dy = -half_sizey + abs(pmax[1])
		
		  print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
		  print, 'Run:' + string(runs)
		  obj_cube[*,*,ii] = shift_sub(obj_cube[*,*,ii], -dx, -dy); shift_sub shifts frames by non-int values, which we're using to center
		endfor; ii for
		
		;Write our centered FITS file
		writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', obj_cube
		print, 'Centered FITS written for run:' + string(runs)
	endfor; runs = 1,2 for left aperture only
	
	if aperture eq 'right' then for runs=3,4 do begin
		; Do this for runs eq 1 and runs eq 3
		if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
		
		output_folder = cube_folder + '/processed_right/'
		
		print, 'Run: ' + string(runs) + 'Output(+dith) Folder: ' + output_folder + dither_folder
		print, 'Reading in a sky_sub cube post-split...'
		obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
		print, 'Post-split sky_sub cube read! Restoring post-split flags and angles...'
		restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
		print, 'Post-split flags and angles read! Changing some stuff in the cube...'
		obj_cube[where(finite(obj_cube) eq 0)] = 0.
		print, 'Some stuff in the cube has been changed! Creating the obj_cube copy...'
		
		;initialize a new cube to modify with the same data
		;obj_cube2 = obj_cube
		;print, 'Obj_cube copy created! Blocking out what needs blocked...'
		
		;we don't need all of the frames, so select some at random
		random_inds=round((n_elements(angles))*randomu(Seed,100))

		obj_cube2=obj_cube[*,*,random_inds]
		angles2=angles[random_inds]


		;first find the peak of the median and shift to there
		medarr,obj_cube2,medframe

		;if jjj eq 1 then writefits,strcompress(reduced_dir+name+'_'+side+'_cube1_align_median_PSF.fits',/rem), medframe
		medframe=medframe-smooth(medframe,25,/nan)
		for xx=0.,2.0*half_cropped_sz-1 do begin
			for yy=0.,2.0*half_cropped_sz-1 do begin
				xd=abs(xx-half_cropped_sz)
				yd=abs(yy-half_cropped_sz)
				if sqrt(xd*xd + yd*yd) gt 50 then medframe[xx,yy]=0.
			endfor
		endfor
		
		writefits,'~/Desktop/centest.fits',medframe
		m=max(medframe,i,/nan)
		ix=i mod (2.0*half_cropped_sz)
		iy=i/(2.0*half_cropped_sz)
		print, 'Center guess = ',ix,iy
		cntrd,medframe,ix,iy,ix,iy,7
		print,'Found center at: ',ix,iy

		dx=float(half_cropped_sz)-ix-0.5
		dy=float(half_cropped_sz)-iy-0.5
		for ii=0,(size(obj_cube))[3]-1 do obj_cube[*,*,ii]=shift_sub(obj_cube[*,*,ii],dx,dy)
;		
;		print, size(obj_cube2)
;		; Divide into oblivion if we're blocking a side out
;		if do_block_right then obj_cube2[fix(0.65*2*half_cropped_sz):fix(2*half_cropped_sz-1), *, *] *= 0.000000000001
;		if do_block_left then obj_cube2[0:fix(0.35*2*half_cropped_sz), *, *] *= 0.000000000001
;		if do_block_bottom then obj_cube2[*, 0:fix(0.35*2*half_cropped_sz), *] *= 0.000000000001
;		if do_block_top then obj_cube2[*,fix(0.65*2*half_cropped_sz):fix(2*half_cropped_sz-1), *] *= 0.000000000001
;		print, 'Requested blocking complete! Changing some stuff in the second cube...'
;		
;		obj_cube2[where(finite(obj_cube) eq 0)] = 0.
;		
;		obj_cube2[where(obj_cube2 lt 0)] *= 0.000000000000001 ;reduce negative values so CC doesn't get confused on subtracted PSF.
;		  print, 'Some stuff in the second cube changed! Cen_filtering if chosen...'
;		half_sizex = half_cropped_sz
;		half_sizey = half_cropped_sz
;		
;		n_frames = (size(obj_cube))[3]
;		if do_cen_filter then begin
;		  for ii=0, n_frames - 1 do begin
;			 print, 'Cen_filter on frame:' + string(ii)
;			 print, 'Run:' + string(runs)
;			 obj_cube2[*,*,ii] -= smooth(obj_cube2[*,*,ii], 17., /NAN)
;			 obj_cube2[*,*,ii] = smooth(obj_cube2[*,*,ii], 5., /NAN)
;		  endfor
;		endif; do_cen_filter if
;		
;		print, 'Cen_filtering complete! (If it was selected). Starting first frame shifting...'
;		; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
;		;Commented out for now, let's try without this first.
;		if do_cen_median_sub then medarr, obj_cube2, medframe
;		
;		for xx=0,2.0*half_sizex-1 do for yy=0,2.0*half_sizey-1 do if sqrt((xx-half_sizex)^2.0+(yy-half_sizey)^2.0) gt align_radius_thresh then obj_cube2[xx,yy,*]=0.
 ;  
;		; take a median and get rid of negatives
;		medarr, obj_cube2, cenmedian
;		cenmedian=cenmedian>0 ;truncate negatives
;
	;	cenmedian=gauss2dfit(cenmedian); fit a gaussian to the median
	;	m=max(cenmedian,i); 
	;	ix=i mod (2.0*half_sizex)
	;	iy=i/(2.0*half_sizey)
	;	cenmedian=shift_sub(cenmedian,half_sizex-ix,half_sizey-iy); shift to center
	;	
	;	for ii=0, n_frames-1 do begin
	;	  ;Commented out for now, let's try without this first.
	;	  if do_cen_median_sub then obj_cube2[*,*,ii] -= medframe
	;	  
	;	  ; Changed from 130
	;	  corr = crosscorr(cenmedian, obj_cube2[*,*,ii], pmax, /cut)
	;	  if pmax[0] ge 0. then dx = half_sizex - pmax[0] else dx = -half_sizex + abs(pmax[0])
	;	  if pmax[1] ge 0. then dy = half_sizey - pmax[1] else dy = -half_sizey + abs(pmax[1])
	;	
	;	  print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
	;	  print, 'Run:' + string(runs)
	;	  obj_cube[*,*,ii] = shift_sub(obj_cube[*,*,ii], -dx, -dy); shift_sub shifts frames by non-int values, which we're using to center
	;	endfor; ii for
	;	
		;Write our centered FITS file
		writefits, output_folder + dither_folder + obj_name +  '_cube_skysub_cen.fits', obj_cube
		print, 'Centered FITS written for run:' + string(runs)
	endfor; runs = 3,4 for right aperture only
	
endif; second stripe
	
print, 'Finished!'

end