pro split, obj, output_folder
;'HII1348', '~/OneDrive/Research/HII1348/testing/' for current testing
COMPILE_OPT IDL2

for runs=1,4 do begin; I think this is the first time that the runs matter...
; Setup (we don't need to read in the obj_cube or angles/oldangles multiple times)
print, 'Starting run: ' + string(runs) + ' , reading fits and retoring flags and angles'
obj_cube=readfits(output_folder+obj+'_cube_skysub.fits')
restore,filename=output_folder+obj+'_parang.sav'
oldangles=temporary(angles)
frames=(size(obj_cube))[3]

   ; Runs 1 and 3
   if runs mod 2 then dither_folder = '/dith1/' else dither_folder = '/dith2/'
   ; Runs 1 and 2
   if runs lt 3 then begin
      side_folder = output_folder + '/processed_left/'
   endif else begin
      side_folder = output_folder + '/processed_right/'
   endelse

print, "Everything is loaded, let's split!"
;fpnod=300./coadd

;I think that this might just need to be manually set depending on which side your first frame is on...
side=0;originally 0, let's try 1?

new_cube=list()
angles=list()


; First frame
print, 'Working on frame: 1 / ' + string(frames)
newframe=obj_cube[*,*,0]
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


print, 'Converting new cube ' + string(runs) + ' to array...'
new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
print, 'New cube converted to array! Cropping...'

if dither_folder eq '/dith1/' then new_cube=new_cube[547-250:547+249,*,*] else new_cube=new_cube[1475-250:1475+249,*,*]
if side_folder eq output_folder + '/processed_left/' then new_cube=new_cube[*,701-250:701+249,*] else new_cube=new_cube[*,259-250:259+249,*] 

print, 'Cropped new cube! Writing FITS...'

writefits, side_folder+dither_folder+obj+'_cube_skysub.fits', new_cube

print, 'FITS written! Converting angles to array...'
angles = angles.toArray(/TRANSPOSE, /NO_COPY)
print, 'Angles converted to array! Saving angles...'

save,filename=side_folder+dither_folder+obj+'_parang.sav',angles
print, 'Angles saved! Finished with run: ' + string(runs)

endfor; runs for
print, 'Done.'
end