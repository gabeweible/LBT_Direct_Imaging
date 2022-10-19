pro create_cube, object_name, raw_fits_path, start_frame, coadd, output_path
; Current run:
; 'HII1348', '~/OneDrive/Research/HII1348/HII1348/raw', 0, 5., '~/OneDrive/Research/HII1348/testing4/'
compile_opt idl2; 32-bit Integers and only square brackets for array indexing
;+
; NAME:
;       CREATE_CUBE
; PURPOSE:
;       Write IDL array frames to a datacube, then write the cube to storage
; CALLING SEQUENCE:
;       CREATE_CUBE, object_name, raw_fits_path, start_frame, coadd, output_path
; INPUTS:
;       object_name = String object name to look for in FITS header
; 
;       raw_fits_path = String containing the path to the raw FITS files
;                   to be added to the cube
;       
;       START_FRAME = Integer frame index to start at in found FITS files
;
;       COADD = Integer number of frames to read before adding to the cube,
;               if coadd is set to 1, then all frames will be added
;
;       output_path = String folder to write the data cube FITS into
;       
; RESTRICTIONS:
;       (1) Make sure that object_name, raw_fits_path, and output_path are string literals
;       and that start_frame and coadd are integers.
;
; EXAMPLE:
;       Write a data cube from FITS files of Vega in ~/Desktop/images starting at
;       frame index 10, with a coadd of 5, outputting to ~/Desktop/Cube:
;       
;       CREATE_CUBE, 'Vega', '~/Desktop/images', 10, 5, '~/Desktop/Cube'
;
; EXTERNAL PROCEDURES USED:
;       READFITS, FXPAR, WRITEFITS
;
; MODIFICATION HISTORY:
;       WRITTEN, Kevin Wagner ????
;       Moved into this procedure and edited by Gabriel Weible 2021
;       Gabriel Weible 2022 changed coadds to be a mean instead of a sum
;-

starttime=systime(/JULIAN)
; Search for the raw images in the specified data path and give some output to the user
print, 'Searching for FITS files in', raw_fits_path, '...'
files = FILE_SEARCH(raw_fits_path, '*.fits', COUNT=filecount)
print, 'Found ', filecount, ' FITS files!'

; Use lists to take advantage of faster appending over arrays (we'll convert back to arrays later)
print, 'Initializing...'
obj_cube = list()
angles = list()
dits = list()
flags = list()

; Initialize an empty coadd frame and our angle for the coadd as 0
coadd_frame = readfits(files[start_frame], head)
replicate_inplace, coadd_frame, 0.
coadd_angle = 0.

; k counts what image we're on, which is used to know when to add the coadded frame to the cube
k = 1
; Start looping through each image
print, 'Initialized successfully, beginning cube creation loop...'
for ii = start_frame, filecount-1 do begin
   print, 'File index', ii, '/', filecount-1
   frame = readfits(files[ii], head)
   angle=fxpar(head, 'LBT_PARA')
   dit = fxpar(head,'ITIME')
   flag = fxpar(head,'FLAG')
         ;perform CDS
         frame = reform(frame[*,*,1] - frame[*,*,0])
         ; Some frames have the wrong dimensions and need cropped
         if (size(frame))[2] gt 2000 then frame = frame[*, 512:1535]; 1535 = 2048 - 512 - 1?
         ; Add to our current coadd frame as a sum, and the same with the current coadd angle
         coadd_frame += frame
         coadd_angle += angle
         if k mod coadd eq 0 then begin
            dits.Add, dit
            flags.Add, flag
            ;Add the mean frame to the cube and reset the frame, and add the mean angle and reset the angle
            coadd_frame *= 1. / coadd
            obj_cube.Add, [[coadd_frame]]
            replicate_inplace, coadd_frame, 0.
            coadd_angle *= 1. / coadd
            angles.Add, coadd_angle
            coadd_angle = 0.
         endif
   print, 'Number of Frames in Cube: ', n_elements(obj_cube)
   k += 1
endfor; ii = start_frame for loop

; Writefits takes an array as input, so we'll need to convert our list over to an array
print, 'Converting to array...'
obj_cube = obj_cube.toArray(/TRANSPOSE, /NO_COPY)
print, 'Writing cube FITS...'

;Write the cube
writefits, strcompress(output_path+object_name+'_cube.fits',/rem), obj_cube
print, 'Cube FITS written! Writing save file...'
angles = angles.toArray(/NO_COPY)
flags = flags.toArray(/NO_COPY)
save, filename=strcompress(output_path+object_name+'_parang.sav',/rem), angles, flags
print, 'Save file written!'
print, 'Completed cube creation in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'

end