pro photo_test
compile_opt IDL2
newline = string(10B)


   output_path = '/Users/gabeweible/OneDrive/research/HII1348/macbook_25/'
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
   half_percent = 15; +- this percent in contrast to grid search
   half_con = half_percent / 100.; Convert to a decimal

   ; Companion centroid guess [x,y] indices (start at 0)
   guess = [277, 353]
   fwhm = 8.72059 ; px ``width'' in reduce_lbti_HII1348.pro
   ; nx x ny grid aroudn centroid result
   nx = 5 & ny = 5
   half_width = 2; px
   
   ; Read in image
   og_image=readfits(strcompress(output_path+'combined/'+obj +'ct_'+string(ct)+'filt_'+$
      string(filter)+'_neg_inj_'+string(0)+'_total_adi.fits',/rem))
   
   ; Do centroiding stuff
   cntrd, og_image, guess[0], guess[1], XCEN, YCEN, fwhm
   print, 'CNTRD xcen, ycen:', XCEN, YCEN
   x_avg = XCEN & y_avg = YCEN

   gcntrd, og_image, guess[0], guess[1], XCEN, YCEN, fwhm
   print, 'GCNTRD xcen, ycen:', XCEN, YCEN

   x_avg += XCEN & x_avg *= 0.5 & y_avg += YCEN & y_avg *= 0.5
   print, 'Mean xcen, ycen:', x_avg, y_avg, newline
   ; This will get us the best place to inject our planet at
   print, 'Starting loop over xx, yy around mean xcen, ycen (nx x ny px square)...'
   
   ; Define lower bounds, upper bounds, and step sizes for our nested loops (grid search)
   x_i = x_avg-half_width & x_f = x_avg+half_width & x_step = ((x_avg+half_width)-(x_avg-half_width))/(nx-1)
   y_i = y_avg-half_width & y_f = y_avg+half_width & y_step = ((y_avg+half_width)-(y_avg-half_width))/(ny-1)
   c_i = c_guess*(1.0-half_con) & c_f = c_guess*(1.0+half_con) & c_step = (c_guess*(1.0+half_con)-c_guess*(1.0-half_con))/(n_contrasts-1)
   
   ; Create arrays to loop through
   x_loop = [x_i : x_f : x_step]
   y_loop = [y_i : y_f : y_step]
   c_loop = [c_i : c_f : c_step]
   
   trial=0
   foreach xx, x_loop do begin; Loop over x
      foreach yy, y_loop do begin; Loop over y

         ; Something isn't working with injecting planets directly, so for now I'm using r and theta
         planet_r = pxscale * SQRT(((xx-250.)^2.)+((yy-250.)^2.))
         planet_theta = ATAN((250.-yy)/(250.-xx)); radians (ccw from top at 0 rad)

         foreach contrast, c_loop do begin; Loop over (negative) contrasts
            print, trial
            trial += 1
         endforeach
      endforeach
   endforeach
   
   help, x_loop
   help, y_loop
   help, c_loop
   
   print, 'x_loop:', x_loop
   print, 'y_loop:', y_loop
   print, 'c_loop:', c_loop
end