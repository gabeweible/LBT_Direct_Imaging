; for rho = ...; inject negative planet to try and cancel out the secondary psf, get the stdev of the
; square around the negative injection to be as low as possible which means that the residuals there
; are as flat as possible, as if the seondary wasn't there (loop over r, theta, and contrast to
; refine it)
;
; aper does aperture photometry, cntrd gets the "center of light" for the sources in the apertures

; I probably need to go full grid-search and loop over negative contrasts at all positions as well.
; I think it makes sense to pass the xx and yy values straight into the planet injection as well,
; as something could go slightly wrong with either conversion to and from radii and angles, but
; I've been having problems with this so for right now the injection is passed through a radius and
; an angle.
; The centroiding should get me within a pixel of the true best negative injection center

pro loop_test, coadd=coadd, type=type
   ; Type is 'ADI' or 'KLIP'

   COMPILE_OPT IDL2
   newline = string(10B)
   start_time = systime(/JULIAN); Get the current time in a Julian date format

   ;------------------------------[ Start User Input ]---------------------------------

   ; Where to find our files and put the results
   output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_'+strcompress(coadd,/r)+'/'
   obj = 'HII1348'

   ; Parameters needed to read in our total_klip or total_adi file
   bin = 3 & bin_type = 'mean' & combine_type = 'nwadi'
   k_klip = 7 & angsep= 1. & anglemax = 360. & nrings = 4.
   n_ang = 2 & do_cen_filter = 1 & filter = 17. & ct = 0.994

   ; Inject_planets parameters
   ; We're only injecting one negative companion at a time
   n_planets = 1
   pxscale = 0.0107 ;arcsec/pixel
   c_guess = -0.0097233772; From sources.txt file (and made negative)
   n_contrasts = 5; Number of contrasts to test at each position

   ; Companion centroid guess [x,y] indices (start at 0)
   guess = [277, 353]
   fwhm = 8.72059 ; ``width'' in reduce_lbti_HII1348.pro
   ; nx x ny grid within centroid pixel
   nx = 5 & ny = 5

   ;------------------------------[ End User Input ]---------------------------------

   ; Initialize arrays for our results
   xxs = [] & yys = [] & cons = [] & devs = [] & means = []

   print, 'Reading in total ' + type +  ' image...'
   if type eq 'KLIP' then begin

      og_image=readfits(strcompress(output_path+'combined/'+obj+'_bin_'+string(sigfig(bin,1))+'_type_'$
         +bin_type+'_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+$
         string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
         string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2))+'neg_inj_'+string(0)+$
         '_total_klip.fits', /rem))

   endif
   if type eq 'ADI' then begin
      og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+'filt_'+$
         string(filter)+'_neg_inj_'+string(0)+'_total_adi.fits',/rem))
   endif
   print, 'Original image read, finding companion centroid...', newline

   ; Calculate the centroid of the companion two ways (at least to test, for now)
   ; Apparently, gcntrd is better, but slower? For now I'll take the average of the two results
   cntrd, og_image, guess[0], guess[1], XCEN, YCEN, fwhm
   print, 'CNTRD xcen, ycen:', XCEN, YCEN
   x_avg = XCEN & y_avg = YCEN

   gcntrd, og_image, guess[0], guess[1], XCEN, YCEN, fwhm
   print, 'GCNTRD xcen, ycen:', XCEN, YCEN

   x_avg += XCEN & x_avg *= 0.5 & y_avg += YCEN & y_avg *= 0.5
   print, 'Mean xcen, ycen:', x_avg, y_avg, newline

   ; Test a 1px x 1px square around the avg centroid position from cntrd and gcntrd
   ; This will get us the best place to inject our planet at
   print, 'Starting loop over xx, yy around mean xcen, ycen (1 px square)...'

   trial = 1 ; Keep track of which run we're on
   for xx = x_avg-0.5, x_avg+0.5, ((x_avg+0.5) - (x_avg-0.5))/(nx-1) do begin
      ;print, 'At xx =', string(xx)
      for yy = y_avg-0.5, y_avg+0.5, ((y_avg+0.5) - (y_avg-0.5))/(ny-1) do begin
         ; print, 'At yy =', string(yy)
         ; Something isn't working with injecting planets directly, so for now I'm using r and theta
         ;planet_r = pxscale * SQRT( ((xx-250.)^2.) + ((yy-250.)^2.) )
         ;planet_theta = ATAN( (250.-yy) / (250.-xx) ); radians (ccw from top at 0 rad)

         for contrast = c_guess*0.95, c_guess*1.05, (c_guess*1.05 - c_guess*0.95) / (n_contrasts-1) do begin
            print, trial

            trial += 1
         endfor; contrast for
      endfor; yy for
   endfor; xx for

   print, 'Loop Test in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'

end; That's all, folks!