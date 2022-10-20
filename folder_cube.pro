pro folder_cube, find_path, save_path
; Searches for all .FITS files in find_path, then combines them into a cube written to save_path

compile_opt IDL2

starttime=systime(/JULIAN)
; Search for the raw images in the specified data path and give some output to the user
print, 'Searching for FITS files in', find_path, '...'
files = FILE_SEARCH(find_path, '*.fits', COUNT=filecount)
print, 'Found ', filecount, ' FITS files!'

; Use lists to take advantage of faster appending over arrays (we'll convert back to arrays later)
print, 'Initializing...'
cube = list()

; k counts what image we're on, which is used to know when to add the coadded frame to the cube
k = 1
; Start looping through each image
print, 'Initialized successfully, beginning cube creation loop...'
for ii = 0, filecount-1 do begin
   print, 'File index', ii, '/', filecount-1
   frame = readfits(files[ii])
   cube.Add, [[frame]]
   print, 'Number of Frames in Cube: ', n_elements(cube)
   k += 1
endfor; ii = start_frame for loop

; Writefits takes an array as input, so we'll need to convert our list over to an array
print, 'Converting to array...'
cube = cube.toArray(/TRANSPOSE, /NO_COPY)
print, 'Writing cube FITS...'

;Write the cube
writefits, strcompress(save_path+'cube.fits',/rem), cube
print, 'Completed cube creation in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'

end