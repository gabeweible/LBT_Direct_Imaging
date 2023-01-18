pro dewarp, dir, obj

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
restore, filename = dir + obj + '_parang.sav'; Restore angles
oldangles = temporary(angles)
frames = (size(obj_cube))[3]; Total # of images in each datacube

; 2 runs because we haven't split the data yet, but do have two mirrors with 
; separate distortion coefficients
for run = 1,2 do begin
   print, 'Starting run: ' + string(run)
   ; different distortion coefficient matrices for each mirror
   if run eq 1 then begin; left mirror (top PSF)
      mirror = 'SX'

      Kx = [[-7.96016109e+00, 9.26584096e-03, -7.93676069e-06, 5.13414639e-10],$
         	[1.02925896e+00, -1.59974177e-05, 9.57769272e-09, -1.14409822e-12],$
         	[-2.30169348e-05, -3.45351550e-09, 1.89621010e-12, 6.72971750e-17],$
         	[7.07041647e-09, 2.20511200e-12, -1.66433082e-15, 2.55463711e-19]]

      Ky = [[-2.26409123e+00, 9.93175401e-01, -6.67169688e-06, 8.08275391e-09],$
         	[-1.38521740e-02, -2.27910031e-05, 4.72367613e-09, -1.64738716e-12],$
         	[8.17060299e-06, 1.35240460e-08, -5.52374318e-12, 2.14966954e-15],$
         	[9.25982725e-10, -2.50607186e-12, 2.58675626e-15, -9.82463036e-19]]
   endif
   
   if run eq 2 then begin; right mirror (bottom PSF)
      mirror = 'DX'

      Kx = [[-1.13034544e+01, 1.45852226e-02, -1.13372175e-05, 1.32602063e-09],$
 	  			[1.03220943e+00, -1.93058980e-05, 1.55798844e-08, -3.86115281e-12],$
 	  			[-2.57352199e-05, 2.70371257e-09, -6.62650011e-12, 3.25017078e-15],$
 	  			[8.02004325e-09, -5.59037685e-13, 1.60256679e-15, -8.18749145e-19]]

      Ky = [[-9.37023860e-01, 9.89973161e-01, -4.32284634e-06, 7.08000564e-09],$
 	  			[-1.59438718e-02, -1.95099497e-05, 1.55801035e-09, 2.13743170e-13],$
 	  			[9.34251811e-06, 1.26473736e-08, -3.71968684e-12, 5.88384784e-16],$
 	  			[3.41071678e-10, -1.82060569e-12, 1.59690189e-15, -3.14810895e-19]]
   endif

   print,"Everything's loaded in - let's dewarp!"
   
   ; Initialize lists for our cube and angles
   new_cube = list()
   warped_cube = list()
   angles = list()
   
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
     angles.Add, oldangles[ii]
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

