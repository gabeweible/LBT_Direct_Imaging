pro dewarp, dir, obj, Kx_sx, Ky_sx, Kx_dx, Ky_dx

; PURPOSE: dewarps LBTI LMIRCAM images using coefficients Kx_dx, Ky_dx, Kx_sx,
; and Ky_sx (update, if needed) for two-mirror, non-interferometric observations,
; resulting in two new cubes being written, one that was dewarped with the sx 
; mirror coefficients, and the other with the dx mirror coefficients, which then
; will need to be differentiated in split.pro.
;
; ARGUMENTS:
; dir: string directory to find the sky-subtraced cube to read in and for 
; writing our dewarped cubes to
;
; obj: string target object name (probably just the host star if directly imaging
; an exoplanet or exoplanets
; (or brown dwarfs! :))
;
; RETURNS: None, but will write two dewarped datacubes to dir (one for the sx
; mirror and one for the dx mirror)
;
; AUTHOR: Gabe Weible, Jan 2023, Undisclosed location

COMPILE_OPT IDL2; Strictarr, 32 bit integers by default
newline = string(10B)

; dir =  '/Users/gabeweible/OneDrive/research/HII1348/macbook_25/'
; obj = 'HII1348'

; Setup
; Read in our sky-subtracted cube to dewarp
obj_cube = readfits(dir + obj + '_cube_skysub.fits')
frames = (size(obj_cube))[3]; Total # of images in each datacube

; 2 runs because we haven't split the data yet, but do have two mirrors with 
; separate distortion coefficients
for run = 1,2 do begin
   print, 'Starting run: ' + string(run)
   ; different distortion coefficient matrices for each mirror
   if run eq 1 then begin; left mirror (top PSF)
      mirror = 'SX'
      Kx = Kx_sx
      Ky = Ky_sx
   endif
   
   if run eq 2 then begin; right mirror (bottom PSF)
      mirror = 'DX'
      Kx = Kx_dx
      Ky = Ky_dx
   endif

   print,"Everything's loaded in - let's dewarp!"
   
   ; Initialize lists for our cube and angles
   new_cube = list()
   warped_cube = list()
   
   for ii=0, frames-1 do begin
     print, 'Working on frame: ' + string(ii + 1) + ' / ' + string(frames)
     print, 'side: ' + mirror + newline
     
     ;Grab a frame from our cube to work on
     warped = obj_cube[*,*,ii]
     padded = replicate(!VALUES.F_NAN, 2048, 2048); full frame
     ; insert warped center stripe into full-frame padding
     ; 512 px (up to index 511) are NaN, 513th px (index 512) is the first with
     ; our actual data and not NaNs
     padded[0, 512] = warped
     
     ; Test before dewarping
     warped_cube.add, [[padded]]
     
     dewarped = POLY_2D(padded, Kx, Ky); dewarp w/coefficient matrices
     
     ; Add to lists
     new_cube.Add, [[dewarped]]
   endfor; ii for

;--------------------------------------------------------------------------------------------------------

   print, 'Converting new cube ' + string(run) + ' to array...'
   new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
   
   ; Test
   warped_cube = warped_cube.toArray(/TRANSPOSE, /NO_COPY)
   
   print, 'New cube converted to array! Writing cube...' + newline
   writefits, dir + obj + '_dewarped_' + mirror + '.fits', new_cube
   
   ; Test
   writefits, dir + obj + '_warped' + '.fits', warped_cube
   
   print, mirror + ' cube written!'

endfor; run for

end

