pro split, obj, output_folder
;'HII1348', '~/OneDrive/Research/HII1348/testing/' for current testing
COMPILE_OPT IDL2
newline = string(10B)

sx_cube = readfits(output_folder + obj + '_dewarped_SX.fits')
dx_cube = readfits(output_folder + obj + '_dewarped_DX.fits')

restore, filename = output_folder + obj + '_parang.sav'
oldangles = temporary(angles)
frames = (size(sx_cube))[3]; Will be the same for both cubes

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
	new_cube = new_cube[547-250:547+249, *, *]; nod on left of frame
endif else BEGIN
	new_cube = new_cube[1475-250:1475+249, *, *]; nod on right of frame
ENDELSE

; Cropping rows (differentiate between the two mirrors)
if side_folder eq output_folder + '/processed_left/' then begin
	new_cube = new_cube[*, 701-250:701+249, *]; top PSF (from SX mirror)
endif else BEGIN
	new_cube = new_cube[*, 259-250:259+249, *]; bottom PSF (from DX mirror)
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
end