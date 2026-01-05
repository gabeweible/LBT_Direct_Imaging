pro split, obj, stripe, output_folder, half_cropped_sz, aperture, skysub_first,$
	do_dewarp, fwhm=fwhm
;'HII1348', '~/OneDrive/Research/HII1348/testing/' for current testing
COMPILE_OPT IDL2
newline = string(10B)

if do_dewarp eq 1 then begin

	if (aperture eq 'left') or (aperture eq 'both') then begin
		sx_cube = readfits_fast(output_folder + obj + '_dewarped_SX.fits')
	endif
	if (aperture eq 'right') or (aperture eq 'both') then begin
		dx_cube = readfits_fast(output_folder + obj + '_dewarped_DX.fits')
	endif
	
endif else begin; do_dewarp if

	if (aperture eq 'left') or (aperture eq 'both') then begin
		if skysub_first eq 0 then print, 'reading skysub cube...' & sx_cube = readfits_fast(output_folder + obj + '_cube_skysub.fits')
		if skysub_first eq 1 then print, 'reading_bad_px_cube...' & sx_cube = readfits_fast(output_folder + obj + '_cube_no_bad_pixels.fits')
	endif
	if (aperture eq 'right') or (aperture eq 'both') then begin
		if skysub_first eq 0 then print, 'reading skysub cube...' & dx_cube = readfits_fast(output_folder + obj + '_cube_skysub.fits')
		if skysub_first eq 1 then print, 'reading_bad_px_cube...' & dx_cube = readfits_fast(output_folder + obj + '_cube_no_bad_pixels.fits')
	endif
	
endelse; not de-warped

if do_dewarp eq 1 then begin
	restore, filename = output_folder + obj + '_parang_dewarp.sav'
endif else begin
	restore, filename = output_folder + obj + '_parang_sky_sub.sav'
endelse
oldangles = temporary(angles)

if aperture eq 'left' then frames = (size(sx_cube))[3]; Will be the same for both cubes
if (aperture eq 'right') or (aperture eq 'both') then frames = (size(dx_cube))[3]

; HII 1348 CENTER 1024 PIXEL ROWS
if (stripe eq 'center1024') and (aperture eq 'both') then begin
	
	for run = 1,4 do begin; Splitting into 4 quadrants (2 mirrors, 2 nods)
	; Setup
	print, 'Starting run: ' + string(run)

	; run 1 and 3
	if run mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
	; run 1 and 2
	if run lt 3 then begin
		; here "left" is the SX mirror, unrelated to the nodding
		side_folder = output_folder + '/processed_left/'
		obj_cube = temporary(sx_cube)
	endif else begin
		; here "right" is the DX mirror, unrelated to the nodding
		side_folder = output_folder + '/processed_right/'
		obj_cube = temporary(dx_cube)
	endelse

	print, "Everything is loaded, let's split!" + newline
	;fpnod=300./coadd

	; I think that this might just need to be manually set depending on which side 
	; the first frame is on (nodding side, not mirror bc we have both all the time)
	side=0;originally 0

	new_cube=list()
	angles=list()

	; First frame
	print, 'Working on frame: 1 / ' + string(frames)
	newframe=obj_cube[*,512:1535,0]

	flag = flags[0]
	next_flag = flags[1]
	if dither_folder eq '/dith1/' and not side then begin
		new_cube.Add, [[newframe]]
		angles.Add, oldangles[0]
	endif
	if dither_folder eq '/dith2/' and side then begin
		new_cube.Add, [[newframe]]
		angles.Add, oldangles[0]
	endif
	if next_flag ne flag then begin
		if side then side=0 else side=1
	endif


	; Middle frames
	for ii=1, frames-2 do begin
	  print, 'Working on frame: ' + string(ii + 1) + ' / ' + string(frames)
	  ;Grab a frame from our cube to work on
	  newframe = obj_cube[*, 512:1535,ii]
  
	  flag = flags[ii]
	  next_flag = flags[ii+1]
	  if dither_folder eq '/dith1/' and not side then begin
			new_cube.Add, [[newframe]]
			angles.Add, oldangles[ii]
	  endif
	  if dither_folder eq '/dith2/' and side then begin
			new_cube.Add, [[newframe]]
			angles.Add, oldangles[ii]
	  endif
	  if next_flag ne flag then begin
		  if side then side=0 else side=1
	  endif
	endfor


	; Last frame
	print, 'Working on frame: ' + string(frames) + ' / ' + string(frames)
	;Grab a frame from our cube to work on
	newframe = obj_cube[*,512:1535,frames-1]

	flag = flags[frames-1]
	if dither_folder eq '/dith1/' and not side then begin
		;print, 'success'
		new_cube.Add, [[newframe]]
		angles.Add, oldangles[frames-1]
	endif; dith 
	if dither_folder eq '/dith2/' and side then begin
		new_cube.Add, [[newframe]]
		angles.Add, oldangles[frames-1]
	endif ;dith2 side = 1 if


	print, 'Converting new cube ' + string(run) + ' to array...'
	new_cube = (temporary(new_cube)).toArray(/TRANSPOSE, /NO_COPY)
	print, 'New cube converted to array! Cropping...'

	; Cropping columns (differentiate between the 2 nods)
	if dither_folder eq '/dith1/' then BEGIN
		new_cube = (temporary(new_cube))[547-half_cropped_sz:547+half_cropped_sz-1, *, *]; nod on left of frame
	endif else BEGIN
		new_cube = (temporary(new_cube))[1475-half_cropped_sz:1475+half_cropped_sz-1, *, *]; nod on right of frame
	ENDELSE

	; Cropping rows (differentiate between the two mirrors)
	if side_folder eq output_folder + '/processed_left/' then begin
		new_cube = (temporary(new_cube))[*, 701-half_cropped_sz:701+half_cropped_sz-1, *]; top PSF (from SX mirror)
	endif else BEGIN
		new_cube = (temporary(new_cube))[*, 259-half_cropped_sz:259+half_cropped_sz-1, *]; bottom PSF (from DX mirror)
	endelse

	print, 'Cropped new cube! Writing FITS...'

	writefits, side_folder+dither_folder+obj+'_cube_skysub.fits', new_cube

	print, 'FITS written! Converting angles to array...'
	angles = (temporary(angles)).toArray(/TRANSPOSE, /NO_COPY)
	print, 'Angles converted to array! Saving angles...'

	save,filename=side_folder+dither_folder+obj+'_parang.sav',angles
	print, 'Angles saved! Finished with run: ' + string(run)

	endfor; run for
	print, 'Done.'

; ALCOR SECOND STRIPE OF 512
endif else if stripe eq 'second512' then begin

	if aperture eq 'left' then begin
		; THESE ARE FROM THE DEWARPED (PADDED) IMAGES OF ALCOR!!!!
		; Also, remember that DS9 starts at (1,1) in the bottom-left corner of a frame
		;(x1, y1) for the left-nod SX PSF, (x2, y2) for the right-nod SX PSF
		SX_PSF_X1 = 440; nod on left of frame
		SX_PSF_X2 = 1364; nod on right of frame
		SX_PSF_Y1 = 1278; px above the bottom of the padded 2048x2048 frame
		SX_PSF_Y2 = 1252
	
		for run = 1,2 do begin; Splitting into 2 sections (1 mirrors, 2 nods)
			; Setup
			print, 'Starting run: ' + string(run)
		
			; run 1 and 3
			if run mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
			
			; here "left" is the SX mirror, unrelated to the nodding
			side_folder = output_folder + '/processed_left/'
			obj_cube = temporary(sx_cube)
		
			print, "Everything is loaded, let's split!" + newline
			;fpnod=300./coadd
		
			; I think that this might just need to be manually set depending on which side 
			; the first frame is on (nodding side, not mirror bc we have both all the time)
			side=0;originally 0; CHANGE THIS FOR ALCOR?? COMPARE WHICH DITHER IS FIRST
			; TO WHAT WAS FIRST FOR HII 1348!!!!!!!
		
			new_cube=list()
			angles=list()
		
		
			; First frame
			print, 'Working on frame: 1 / ' + string(frames)
			newframe=obj_cube[*,*,0]
		
			; crop back to second stripe (right from left, up from bottom)
			;newframe = newframe[*, 1024:1535]
		
			flag = flags[0]
			next_flag = flags[1]
			if dither_folder eq '/dith1/' and not side then begin
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[0]
			endif
			if dither_folder eq '/dith2/' and side then begin
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[0]
			endif
			if next_flag ne flag then begin
				if side then side=0 else side=1
			endif
		
		
			; Middle frames
			for ii=1, frames-2 do begin
			  print, 'Working on frame: ' + string(ii + 1) + ' / ' + string(frames)
			  ;Grab a frame from our cube to work on
			  newframe = obj_cube[*,*,ii]
		  
			 ; crop back to second stripe
			 ;newframe = newframe[*, 1024:1535]
		  
			  flag = flags[ii]
			  next_flag = flags[ii+1]
			  if dither_folder eq '/dith1/' and not side then begin
					new_cube.Add, [[newframe]]
					angles.Add, oldangles[ii]
			  endif
			  if dither_folder eq '/dith2/' and side then begin
					new_cube.Add, [[newframe]]
					angles.Add, oldangles[ii]
			  endif
			  if next_flag ne flag then begin
				  if side then side=0 else side=1
			  endif
			endfor
		
		
			; Last frame
			print, 'Working on frame: ' + string(frames) + ' / ' + string(frames)
			;Grab a frame from our cube to work on
			newframe = obj_cube[*,*,frames-1]
		
			; crop back to second stripe
			;newframe = newframe[*, 1024:1535]
		
			flag = flags[frames-1]
			if dither_folder eq '/dith1/' and not side then begin
				;print, 'success'
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[frames-1]
			endif; dith 
			if dither_folder eq '/dith2/' and side then begin
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[frames-1]
			endif ;dith2 side = 1 if
		
		
			print, 'Converting new cube ' + string(run) + ' to array...'
			new_cube = (temporary(new_cube)).toArray(/TRANSPOSE, /NO_COPY)
			print, 'New cube converted to array! Cropping...'
		
			; Left PSF
			if dither_folder eq '/dith1/' then BEGIN
				; Cropping columns (differentiate between the 2 nods)
				new_cube = (temporary(new_cube))[SX_PSF_X1-half_cropped_sz:SX_PSF_X1+half_cropped_sz-1, *, *]; nod on left of frame
				; Cropping rows
				new_cube = (temporary(new_cube))[*, SX_PSF_Y1-half_cropped_sz:SX_PSF_Y1+half_cropped_sz-1, *]; top PSF (from SX mirror)
				
			endif else BEGIN
				; Cropping columns (differentiate between the 2 nods)
				new_cube = (temporary(new_cube))[SX_PSF_X2-half_cropped_sz:SX_PSF_X2+half_cropped_sz-1, *, *]; nod on right of frame
				;Cropping rows
				new_cube = (temporary(new_cube))[*, SX_PSF_Y2-half_cropped_sz:SX_PSF_Y2+half_cropped_sz-1, *]
			ENDELSE
			print, 'Cropped new cube! Writing FITS...'
		
			writefits, side_folder+dither_folder+obj+'_cube_skysub.fits', new_cube
		
			print, 'FITS written! Converting angles to array...'
			angles = (temporary(angles)).toArray(/TRANSPOSE, /NO_COPY)
			print, 'Angles converted to array! Saving angles...'
		
			save,filename=side_folder+dither_folder+obj+'_parang.sav',angles
			print, 'Angles saved! Finished with run: ' + string(run)
		endfor; run1 and run2 for
		print, 'Done.'
	endif; left
	
	if aperture eq 'right' then begin
		; THESE ARE FROM THE DEWARPED (PADDED) IMAGES OF ALCOR!!!!
		; Also, remember that DS9 starts at (1,1) in the bottom-left corner of a frame
		;(x1, y1) for the left-nod SX PSF, (x2, y2) for the right-nod SX PSF
		if do_dewarp eq 1 then begin
			DX_PSF_X1 = 447; nod on left of frame
			DX_PSF_X2 = 1370; nod on right of frame
			DX_PSF_Y1 = 1281; px above the bottom of the padded 2048x2048 frame
			DX_PSF_Y2 = 1260
		endif else begin; dewarping if
			DX_PSF_X1 = 435; nod on left of frame
			DX_PSF_X2 = 1353; nod on right of frame
			DX_PSF_Y1 = 259; px above the bottom of the padded 2048x2048 frame
			DX_PSF_Y2 = 230
		endelse; not dewarping else
	
		for run = 3,4 do begin; Splitting into 2 sections (1 mirrors, 2 nods)
			; Setup
			print, 'Starting run: ' + string(run)
		
			; run 1 and 3
			if run mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
			
			; here "right" is the DX mirror, unrelated to the nodding
			side_folder = output_folder + '/processed_right/'
			obj_cube = dx_cube
		
			print, "Everything is loaded, let's split!" + newline
		
			; I think that this might just need to be manually set depending on which side 
			; the first frame is on (nodding side, not mirror bc we have both all the time)
			side=0;originally 0; CHANGE THIS FOR ALCOR?? COMPARE WHICH DITHER IS FIRST
			; TO WHAT WAS FIRST FOR HII 1348!!!!!!!
		
			new_cube=list()
			angles=list()
		
			; First frame
			print, 'Working on frame: 1 / ' + string(frames)
			newframe=obj_cube[*,*,0]
		
			; crop back to second stripe (right from left, up from bottom)
			;newframe = newframe[*, 1024:1535]
		
			flag = flags[0]
			next_flag = flags[1]
			if dither_folder eq '/dith1/' and not side then begin
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[0]
			endif
			if dither_folder eq '/dith2/' and side then begin
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[0]
			endif
			if next_flag ne flag then begin
				if side then side=0 else side=1
			endif
		
		
			; Middle frames
			for ii=1, frames-2 do begin
			  print, 'Working on frame: ' + string(ii + 1) + ' / ' + string(frames)
			  ;Grab a frame from our cube to work on
			  newframe = obj_cube[*,*,ii]
		  
			  flag = flags[ii]
			  next_flag = flags[ii+1]
			  if dither_folder eq '/dith1/' and not side then begin
					new_cube.Add, [[newframe]]
					angles.Add, oldangles[ii]
			  endif
			  if dither_folder eq '/dith2/' and side then begin
					new_cube.Add, [[newframe]]
					angles.Add, oldangles[ii]
			  endif
			  if next_flag ne flag then begin
				  if side then side=0 else side=1
			  endif
			endfor
		
		
			; Last frame
			print, 'Working on frame: ' + string(frames) + ' / ' + string(frames)
			;Grab a frame from our cube to work on
			newframe = obj_cube[*,*,frames-1]
		
		
			flag = flags[frames-1]
			if dither_folder eq '/dith1/' and not side then begin
				;print, 'success'
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[frames-1]
			endif; dith 
			if dither_folder eq '/dith2/' and side then begin
				new_cube.Add, [[newframe]]
				angles.Add, oldangles[frames-1]
			endif ;dith2 side = 1 if
		
		
			print, 'Converting new cube ' + string(run) + ' to array...'
			new_cube = (temporary(new_cube)).toArray(/TRANSPOSE, /NO_COPY)
			print, 'New cube converted to array! Cropping...'
		
			; Right PSF
			if dither_folder eq '/dith1/' then BEGIN
				; Cropping columns (differentiate between the 2 nods)
				new_cube = (temporary(new_cube))[DX_PSF_X1-half_cropped_sz:DX_PSF_X1+half_cropped_sz-1, *, *]; nod on left of frame
				; Cropping rows
				new_cube = new_cube[*, DX_PSF_Y1-half_cropped_sz:DX_PSF_Y1+half_cropped_sz-1, *]; top PSF (from SX mirror)
				
			endif else BEGIN
				; Cropping columns (differentiate between the 2 nods)
				new_cube = (temporary(new_cube))[DX_PSF_X2-half_cropped_sz:DX_PSF_X2+half_cropped_sz-1, *, *]; nod on right of frame
				;Cropping rows
				new_cube = (temporary(new_cube))[*, DX_PSF_Y2-half_cropped_sz:DX_PSF_Y2+half_cropped_sz-1, *]
			ENDELSE
		
			print, 'Cropped new cube! Writing FITS...'
		
			writefits, side_folder+dither_folder+obj+'_cube_skysub.fits', new_cube
		
			print, 'FITS written! Converting angles to array...'
			angles = (temporary(angles)).toArray(/TRANSPOSE, /NO_COPY)
			print, 'Angles converted to array! Saving angles...'
		
			save,filename=side_folder+dither_folder+obj+'_parang.sav',angles
			print, 'Angles saved! Finished with run: ' + string(run)
	
		endfor; run for
		print, 'Done.'
	endif; right
	
endif

end