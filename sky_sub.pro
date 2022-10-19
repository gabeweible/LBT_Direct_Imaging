pro sky_sub, obj_name, coadd, output_folder, cube_folder=cube_folder, fs_start=fs_start, fs_end=fs_end
; 'HII1348', 5., '~/OneDrive/Research/HII1348/testing/' for current run (same output and cube folders)
compile_opt idl2

newline = string(10B); Newline character to separate out print statements

if keyword_set(cube_folder) then cube_folder=cube_folder else cube_folder=output_folder

obj_cube = readfits(cube_folder + obj_name + '_cube_no_bad_pixels.fits')
      ; restores angles and flags from the combining section (why do we need to restore them?) Just if we want to first combine into a cube and then do sky sub later?
restore, cube_folder + obj_name + '_parang.sav'
;frames in each nod = 300/bin, e.g. bin=50 => 6 frames in each nod
frames_per_nod = 200. / coadd
      
n_frames = (size(obj_cube))[3]
print, 'Number of frames = ', n_frames
print, 'Number of nods = ', n_frames / frames_per_nod
  
firstsky = obj_cube[*,*,fs_start:fs_end]
firstsky = mean(firstsky, dim=3)
writefits,'~/Desktop/test_firstsky.fits', firstsky
skys = list()
  
sky = temporary(firstsky)

;First frame
flag_i = flags[0]
flag_next_i = flags[1]
print, 'Frame index = ', 0, '/', n_frames-1
print, 'This, next flag = ', flag_i, flag_next_i
;if ii ge 40 then frames_per_nod = 400. / coadd
frame = obj_cube[*,*,0]
skys.Add, [[frame]]
frame -= sky
; Assign our frame as the new cube if we're on the first frame, or add it to the array if we're on any subsequent frame
new_cube = list(frame)


; Middle frames
for ii=1, n_frames-2 do begin
  flag_i = flags[ii]
  flag_next_i = flags[ii+1]
  print, 'Frame index = ', ii, '/', n_frames-1
  print, 'This, next flag = ', flag_i, flag_next_i, newline
  ;if ii ge 40 then frames_per_nod = 400. / coadd
  frame = obj_cube[*,*,ii]
  skys.Add, [[frame]]
  if flag_next_i ne flag_i then begin
   sky = median(skys.toArray(/TRANSPOSE, /NO_COPY), dim=3)
   skys = list()
  endif
   ; This puts the subtraction in sky subtraction!
  frame -= sky
   ; Assign our frame as the new cube if we're on the first frame, or add it to the array if we're on any subsequent frame
  new_cube.Add, [[frame]]
endfor; sky subtraction for

; Last frame
flag_i = flags[n_frames-1]
print, 'Frame index = ', n_frames-1, '/', n_frames-1
print, 'This, next flag = ', flag_i, flag_next_i
;if ii ge 40 then frames_per_nod = 400. / coadd
frame = obj_cube[*,*,n_frames-1]
skys.Add, [[frame]]
; This puts the subtraction in sky subtraction!
frame -= sky
; Assign our frame as the new cube if we're on the first frame, or add it to the array if we're on any subsequent frame
new_cube.Add, [[frame]]

  

;;Save the new_cube list
;save, filename=strcompress(output_folder+obj_name+'_new_cube_list.sav',/rem), new_cube

;print, 'New cube list saved! Converting to array and reassigning...'
;/TRANSPOSE will give us the correct dimensions, /NO_COPY will just move over the list elements to the array, not copying them

print, 'Converting new cube to array...'
new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
print, 'New cube converted to array! Writing FITS file...'
writefits, output_folder + obj_name + '_cube_skysub.fits', new_cube
print, 'FITS file written!'
print, 'Finished!'

end