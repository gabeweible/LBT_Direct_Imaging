pro split, obj, stripe, output_folder, half_cropped_sz
;'HII1348', '~/OneDrive/Research/HII1348/testing/' for current testing
COMPILE_OPT IDL2
newline = string(10B)

sx_cube = readfits(output_folder + obj + '_dewarped_SX.fits')

restore, filename = output_folder + obj + '_parang.sav'
oldangles = temporary(angles)
frames = (size(sx_cube))[3]; Will be the same for both cubes

; HII 1348 CENTER 1024 PIXEL ROWS
if stripe eq 'center1024' then begin

	; only need the DX cube if we're doing center-stripe double-sided obs (or just DX)
	dx_cube = readfits(output_folder + obj + '_dewarped_DX.fits')
	
	for run = 1,4 do begin; Splitting into 4 quadrants (2 mirrors, 2 nods)
	; Setup
	print, 'Starting run: ' + string(run)

	; run 1 and 3
	if run mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
	; run 1 and 2
	if run lt 3 then begin
		; here "left" is the SX mirror, unrelated to the nodding
		side_folder = output_folder + '/processed_left/'
		obj_cube = sx_cube
	endif else begin
		; here "right" is the DX mirror, unrelated to the nodding
		side_folder = output_folder + '/processed_right/'
		obj_cube = dx_cube
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
	newframe=obj_cube[*,*,0]

	; crop back to center stripe
	newframe = newframe[*, 512:1535]

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
  
	  ; crop back to center stripe
	  newframe = newframe[*, 512:1535]
  
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

	; crop back to center stripe
	newframe = newframe[*, 512:1535]

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
	new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
	print, 'New cube converted to array! Cropping...'

	; Cropping columns (differentiate between the 2 nods)
	if dither_folder eq '/dith1/' then BEGIN
		new_cube = new_cube[547-half_cropped_sz:547+half_cropped_sz-1, *, *]; nod on left of frame
	endif else BEGIN
		new_cube = new_cube[1475-half_cropped_sz:1475+half_cropped_sz-1, *, *]; nod on right of frame
	ENDELSE

	; Cropping rows (differentiate between the two mirrors)
	if side_folder eq output_folder + '/processed_left/' then begin
		new_cube = new_cube[*, 701-half_cropped_sz:701+half_cropped_sz-1, *]; top PSF (from SX mirror)
	endif else BEGIN
		new_cube = new_cube[*, 259-half_cropped_sz:259+half_cropped_sz-1, *]; bottom PSF (from DX mirror)
	endelse

	print, 'Cropped new cube! Writing FITS...'

	writefits, side_folder+dither_folder+obj+'_cube_skysub.fits', new_cube

	print, 'FITS written! Converting angles to array...'
	angles = angles.toArray(/TRANSPOSE, /NO_COPY)
	print, 'Angles converted to array! Saving angles...'

	save,filename=side_folder+dither_folder+obj+'_parang.sav',angles
	print, 'Angles saved! Finished with run: ' + string(run)

	endfor; run for
	print, 'Done.'

; ALCOR SECOND STRIPE OF 512
endif else if stripe eq 'second512' then begin

	; THESE ARE FROM THE DEWARPED (PADDED) IMAGES OF ALCOR!!!!
	; Also, remember that DS9 starts at (1,1) in the bottom-left corner of a frame
	;(x1, y1) for the left-nod SX PSF, (x2, y2) for the right-nod SX PSF
	SX_PSF_X1 = 440; nod on left of frame
	SX_PSF_X2 = 1362; nod on right of frame
	SX_PSF_Y1 = 1297; px above the bottom of the padded 2048x2048 frame
	SX_PSF_Y2 = 1270

	for run = 1,2 do begin; Splitting into 2 sections (1 mirrors, 2 nods)
	; Setup
	print, 'Starting run: ' + string(run)

	; run 1 and 3
	if run mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
	
	; here "left" is the SX mirror, unrelated to the nodding
	side_folder = output_folder + '/processed_left/'
	obj_cube = sx_cube

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
	new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
	print, 'New cube converted to array! Cropping...'

	; Left PSF
	if dither_folder eq '/dith1/' then BEGIN
		; Cropping columns (differentiate between the 2 nods)
		new_cube = new_cube[SX_PSF_X1-half_cropped_sz:SX_PSF_X1+half_cropped_sz-1, *, *]; nod on left of frame
		; Cropping rows
		new_cube = new_cube[*, SX_PSF_Y1-half_cropped_sz:SX_PSF_Y1+half_cropped_sz-1, *]; top PSF (from SX mirror)
		
	endif else BEGIN
		; Cropping columns (differentiate between the 2 nods)
		new_cube = new_cube[SX_PSF_X2-half_cropped_sz:SX_PSF_X2+half_cropped_sz-1, *, *]; nod on right of frame
		;Cropping rows
		new_cube = new_cube[*, SX_PSF_Y2-half_cropped_sz:SX_PSF_Y2+half_cropped_sz-1, *]
		ENDELSE

	;TESTING!!!! REMOVE IF WORKING CORRECTLY
	writefits, '~/Desktop/Alcor_crop_test.fits', new_cube
	
	; Cropping rows (differentiate between the two mirrors)
	; approx. avg. height of SX PSF is 269 px above the bottom of the
	; second512 frame. 
	; just remove some white from the bottom. There are only 


	print, 'Cropped new cube! Writing FITS...'

	writefits, side_folder+dither_folder+obj+'_cube_skysub.fits', new_cube

	print, 'FITS written! Converting angles to array...'
	angles = angles.toArray(/TRANSPOSE, /NO_COPY)
	print, 'Angles converted to array! Saving angles...'

	save,filename=side_folder+dither_folder+obj+'_parang.sav',angles
	print, 'Angles saved! Finished with run: ' + string(run)

	endfor; run for
	print, 'Done.'
endif

end