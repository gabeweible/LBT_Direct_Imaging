pro my_HII1348, theta=theta, rho=rho, contrast=contrast
COMPILE_OPT IDL2; Use only square brackets for array indices, not parentheses, and default to 32 bit integers (reduces chance of overflow)
   ;theta, rho, and contrast are keywords (can be passed or not--optional for planet injection)
   start_time = systime(/JULIAN); Get the current time in a Julian date format (number of days since January 1, 4713 BC, with decimals)
   
   ;------------------------------[ Start User Input ]---------------------------------
   
   ;what object name should we look for in the header of our FITS files?
   obj_name = 'HII 1348'
   ; paths for where to find our raw data and for where to put our processed (uncombined) data
   raw_path = '/home/ggdoubleu/OneDrive/Research/HII1348/HII1348/raw'
   left_output_path = '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed_left/'
   right_output_path = '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed_right/'
   
   coadd = 1  ;25;should probably leave this set to 1 to avoid combining jitter, bad frames, etc. (I'll take your word for it and set it to 1!)
   
   ;True/False options: 0 = False, 1 = True (hint: look for the "do_" at the beginning)
   unsat = 0; Required to be 0 for do_combine as written
   do_combine = 1; combine raw data files into a single cube for easier handling
   do_bad_pixels = 0
   do_sky_subtraction = 0; Note: this is required for do_split as written
   do_split = 0 ;splits the dithers into two datasets
   do_center = 0
   do_cen_filter = 0;filters during the centering step (leaves science cube alone)
   do_cen_median_sub = 0
   do_block_right = 0;will cross correlate the left half of the array to not confuse beams
   do_block_left = 0
   do_block_top = 0
   do_block_bottom = 0
   do_clean = 0
   use_injection = 0

   do_block_center = 0 ;WARNING! Will overwrite all centering steps.
   inner_rad = 0.3 ;arcsec
   do_block_outer = 0  ;WARNING! Will overwrite all centering steps.
   outer_rad = 3.5 ;arcsec
   
   ;Number of planets and their radii
   n_planets = 3
   planet_r = [0.6, 0.6, 0.6];findgen(n_planets) * 0.1 + 0.1;[0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

   ; Why such weird multiples of pi being added/multiplied? Is there a purpose to that? This makes 3 evenly spaced planets, for example
   planet_theta = findgen(n_planets, START=1) * 2 * !DPI / n_planets; + !DPI + !DPI/8 - !DPI/2 + !DPI/12 ;+ !DPI
   planet_contrast = fltarr(n_planets) + (2E-5) * 0.5 ;0.5 is for correcting number of coadds
   pixel_scale = 0.0107  ;arcsec/pixel

   ;What are these?
   do_annmode = 0
   annmode_inout = [10, 95]
   ;annmode_inout=[5,25]

   ; Set the necessary values if rho, theta, and or contrast are passed to the procedure for planet injection
   if keyword_set(rho) then begin
     do_inject = 1
     planet_theta = [theta]
     planet_r = [rho]
     use_injection = 1
     do_annmode = 1
     annmode_inout = round([max([0., (planet_r / pixel_scale) - 9.]), (planet_r / pixel_scale) + 12.])
     if annmode_inout[1] gt 16 then annmode_inout = round([max([0.,planet_r / pixel_scale-12.+1.]),planet_r / pixel_scale+12.+2.])
   endif else do_inject = 0

   if keyword_set(contrast) then planet_contrast = [contrast]
   n_planets = n_elements(planet_r)
   
   ;next three parameters used for both rotate and KLIP
   do_destripe = 0  ;will also be used in the precise centering
   filter = 15.

   ;What is nwadi?
   combine_type = 'nwadi' ;'mean','median','nwadi'
   do_rotate = 0
   do_klip = 0
   bin = 4; Is this 2 x 2 = 4 or 4 x 4?
   bin_type = 'mean' ;median or mean
   klip_fraction = 0 ; Does this mean only doing KLIP on a fraction of the frames?
   start_frame = 75
   end_frame = 150
   do_hyper = 0  ;will only process the region around the planet
   fill = 0  ;will fill in the rest of the image with the filtered image
   k_klip = 6 ;max([1,fix(3.*randomu(Seed))])
   ang_sep = 1. ;2.*randomu(Seed)
   max_angle = 360. ;defined in degrees. Defines max angle
   n_rings = 4. ;4;max([1,fix(12.*randomu(Seed))])
   wr = fix(250. / float(n_rings))
   n_angle = 6 ;max([1,fix(9.*randomu(Seed))])

   do_rdi_klip = 0
   rdi_high_pass = 15.
   ref_folder = '/Volumes/Storage/LBT/MWC758/oct_2016_data/processed-ref/'; Do I need this?
   ref_name = 'RY Tau  '
   do_cross_corr_sci = 0 ;cross correlation shift reference cube to the science cube

   if do_rdi_klip then begin
     klip_fraction = 1
     start_frame = 0
     end_frame = 10
     do_hyper = 0  ;will only process the region around the planet
     fill = 1 ;will fill in the rest of the image with the filtered image
     ;k_klip=1;max([1,fix(3.*randomu(Seed))])
     n_rings = 4;max([1,fix(12.*randomu(Seed))])
     wr = fix(140. / float(n_rings))
     n_angle = 4;max([1,fix(9.*randomu(Seed))])
   endif; do_rdi_klip if

   suffix = strcompress(reform('_'+string(k_klip)+'k_'+string(sigfig(ang_sep,2))+'as_'+string(sigfig(max_angle,2))+'am_'+string(sigfig(n_rings,2))+'rings_'+string(wr)+'wr_'+String(n_angle)+'nang_'+string(sigfig(filter,2))+'filter_'+string(sigfig(bin,2))+'bin'), /remove_all)
   ;override suffix defined previously
   ;suffix='_as10_am360_k2_hpf25_n_rings2';strcompress('_combine-'+combine_type,/rem)

   if use_injection then suffix += '_inj'; I don't think we need the strcompress here since we already did it above and '_inj' has no spaces
   ;if do_hyper then suffix=strcompress(reform('_'+string(k_klip)+'k_'+string(sigfig(ang_sep,2))+'as_'+string(sigfig(n_rings,2))+'rings_'+string(wr)+'wr_'+String(n_angle)+'nang_'+string(sigfig(high_pass_width,2))+'filter'+string(sigfig(low_pass_width,2))+'lowfilter_do_hyperer'), /remove_all)
  
   ; Initialize some empty arrays
   ang_seps = []
   max_angles = []
   k_klips = []
   n_angles = []
   bins = []
   filters = []
   SNRs = []
  
   ;for ang_sep = 0.1, 0.1, 0.1 do begin
   ;for max_anglex = 20, 180, 20 do begin
   ;for k_klip = 5, 50, 5 do begin
   ;for n_angle = 2, 6, 2 do begin
   ;for bin = 8,8,1 do begin
   ;for filter = 0., 0., 50. do begin
  
   ;k_klip = 5
   ;n_angle = 6
  
   ;ang_sep = 0.5
   ;max_anglex = 360.
   ;k_klip = 10.
   ;n_angle = 6.
   ;bin = 4
   ;filter = 15.
  
   ;run a loop of 4 different reductions (2 nods * 2 mirrors)
   for runs = 1,4 do begin
      ;raw_path to data
      ;raw_path = '/Volumes/Storage/LBT/MWC758/oct_2016_data/copies_data_files/'
      ;output_folder = '/Volumes/Storage/LBT/MWC758/oct_2016_data/processed-3/'

      ; Runs 1 and 2 are for the left mirror, 3 and 4 are for the right mirror (2 nods/dithers? each)
      if runs eq 1 or runs eq 2 then begin
         output_folder = left_output_path
         side = 'left'
      endif else begin
         output_folder = right_output_path
         side = 'right'
      endelse
    
      ;I don't have these, but probably don't need them?
      ;if runs eq 1 or runs eq 2 then raw_path='/home/ggdoubleu/OneDrive/Research/HII1348/Testing/raw-dewarp-left/'
      ;if runs eq 3 or runs eq 4 then raw_path='/home/ggdoubleu/OneDrive/Research/HII1348/Testing/raw-dewarp-right/'
    
      ;rest of the scripts require the correct dither folder to be used
      if runs eq 1 or runs eq 3 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
      
      ;if output_folder eq '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed-right/' then begin
         ;if dither_folder eq '/dith1/' then first_center=[238,149.]-0.5;[373.5,643.5]-0.5
         ;if dither_folder eq '/dith2/' then first_center=[705.,168.]-0.5;[381.,218.]-0.5
      ;endif
    	
      ;if output_folder eq '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed-left/' then begin   
         ;if dither_folder eq '/dith1/' then first_center=[241,166.]-0.5;[373.5,643.5]-0.5
         ;if dither_folder eq '/dith2/' then first_center=[704.,175.]-0.5;[381.,218.]-0.5
      ;endif

      corr_thresh = 0.99
      if runs ne 1 then corr_thresh = corr_thresh; What does this do? corr_thresh never seems to change
    
      ;do_inject = 0
      ;if do_inject then do_filter = 1	;ensures fresh injection each time
      ;I don't have this, do I need it?
      ;ref_file = '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/unsat/MWC758_PSF_151.085ms_50coadd.fits'
      ;ref file not needed given unsaturated data
    
      ;n_planets=1
      ;planet_r=[0.11];findgen(n_planets)*0.1+0.1;[0.2,0.3,0.4,0.5,0.6,0.7]
      ;planet_theta= (165.+90.+180.)*!DTOR ; findgen(n_planets) * !DPI / (n_planets/2.) + 4. * !DPI / 4. + !DPI/8 - !DPI/2 + !DPI/12 ;+ !DPI
      ;planet_contrast=fltarr(n_planets)
      ;planet_contrast[*]=(10.0^(-6.5/2.5)) * 0.5 ;0.5 is for correcting number of coadds
      ;pixel_scale=0.0107	;arcsec/pixel
    
      ;make a bigger array of fake planets
      ;n_planets = 1
      ;planet_r = [planet_r, planet_r + 0.35]
      ;planet_theta = [planet_theta, planet_theta + !DPI / 8]
      ;planet_contrast = fltarr(n_planets)
      ;planet_contrast[*] = 1.0E-5
    
      ;make a spiral of planets
      ;n_planets = 24
      ;planet_r = findgen(24)*0.25
      ;planet_theta = findgen(24)*!DPI/8
      ;planet_contrast = fltarr(n_planets)
      ;planet_contrast[*] = 1E-5
        
      ;make a weird distribution of planets
      ;n_planets = 48
      ;planet_r = findgen(48) * 0.05
      ;planet_theta = findgen(48) * !DPI / 8 * sin( findgen(48) * !DPI / 4 )
      ;planet_contrast = fltarr(n_planets)
      ;planet_contrast[*] = (0.5E-5); * sin( findgen(48) * !DPI / 8 ) + 1E-5
    
      ;------------------------------[ End User Input ]---------------------------------
      
      ;------------------------------[ Begin Thought Process (Combining Into a Cube) ]---------------------------------
    
      ; Here we combine the individual raw FITS files from the observation into a data cube of frames
      if do_combine then begin
         print, 'Path to raw files:', raw_path
         files = FILE_SEARCH(raw_path, '*.fits', COUNT=filecount)
         print, 'Found ', filecount, ' FITS files.'
         ;hak
         print, 'Assembling data cube...'
         
    	   obj_cube = []
    	   angles = []
    	   k = 0
    	   integ_times = []
    	   flags = []
    	
         start = 4473 ;why do we start here specifically?
         for ii = start, filecount-1 do begin
            frame = readfits(files[ii], head)
    
            ;if side eq 'left' then frame = frame[450:1450, 550:900, *] 
    		    ;if side eq 'right' then frame = frame[450:1450, 150:500, *]
            ;frame = frame[*, 512:1535]
            ;writefits,'~/Desktop/test.fits',frame
    
            ; Do this stuff for the first FITS file only
    		    if ii eq start then begin
    		       framei = frame
    		       framei[*] = 0.
    		       anglei = 0.
    		    endif
    		
    		    ; Get the object name and parallactic angle from the FITS file
    	      object = fxpar(head, 'OBJNAME')
    	      angle = fxpar(head, 'LBT_PARA')
    	
    	      ; make sure that the object name in the FITS file header is a string and matches what we expect
    	      ;PROBLEM: ISA can only use STRING as a keyword in IDL 8.4 and newer, I have 8.3 (replace with typename, should do the same thing in this case)
    	      if typename(object) eq "STRING" and object eq obj_name then begin
               ;if object ne obj_name then begin & print, object & hak & endif
  		         ;print, 'Match found!'
  		         ;if unsat and fxpar(head,'LMIR_FW2') eq 'ND2.0-T1' then begin
  		        
  		         ;Get more info from our FITS file header
    	         if not unsat and fxpar(head,'LMIR_FW2') ne 'ND2.0-T1' then begin; What would if mean if it did equal ND2.0-T1?
    		          ;writefits,'/Volumes/Storage/LBT/test.fits',frame,head
    		          k += 1
    		          integ_time = fxpar(head, 'ITIME')
    		          flag = fxpar(head, 'FLAG')
    		          ;if ii mod 10 eq 0 then begin
    		             
    		          ;endif
    		          ;if k then integ_times = integ_time else 
    		          ;if k then flags = flag else 
    		          ;if total(integ_times-mean(integ_time)) ne 0 then begin print, 'Different integrations found! Reprogram me please to deal with this.' & hak & endif
      	          ;print, integ_time / 151.085
      	       
      	          ; Test to make sure that the integration time for each FITS file makes sense
      	          if integ_time / 906.51 ge 0.95 and  integ_time / 906.51 le 1.05 then begin
                     ;perform CDS (what's that?)
      		           first_frame = frame[*,*,0]
      		           last_frame = frame[*,*,1]
      		           frame = reform(last_frame - first_frame)
      		
      		           ; this line is mysterious to me
        	           if (size(frame))[2] gt 2000 then frame = frame[*, 512:1535]; 1535 = 2048 - 512 - 1?
        	
        	           framei += frame
        	           ;hak
        	
        	           anglei += angle
        	           ; if coadd = 1 as suggested, a mod 1 = 0 for all integer a, so we add all of our frames to the cube. I think that coadd being 25 means we only add every 25th frame to the cube
                     if k mod coadd eq 0 then begin
                        ;framei = framei / float(coadd)
            	          obj_cube = [ [[obj_cube]], [[framei]] ]; add our frame to the cube
            	          integ_times = [integ_times, integ_time]
            	          flags = [flags, flag]
            	          framei[*] = 0.; Reset the frame so that we can start again
            	          anglei = anglei / float(coadd)
            	          angles = [angles, angle]
            	          anglei = 0.
                     endif; k mod coadd eq 0
                    
                  endif ;integ_time if
                 
               endif ;unsat ND filt if
              
            endif	else print, 'No match' ; FITS header name is a string and is what's expected if
            
            cls; clear the IDL terminal, not working right now...
            print, 'Frame index:' + string(ii) + ' /' + string(filecount - 1)
            print, 'Object : ', object
            print, 'Run Number: ' + STRING(runs)
            print, fxpar(head, 'LMIR_FW2')
            print, 'Frame sizes: ', size(frame)
            print, 'Integration time: ', integ_time
            ;flag indices will be off if binning!
            print, 'frame in coadd:', k mod coadd
            print, 'coadd:', coadd
            print, 'Cube size:', size(obj_cube)
           
         endfor ;ii for (for each raw FITS image)
         ; Write the data cube!  
         writefits, strcompress(output_folder + obj_name + '_cube.fits', /rem), obj_cube
         save, filename = strcompress(output_folder + obj_name + '_parang.sav', /rem), angles, flags
        
      endif ;do_combine if
      ;n_frames = (size(obj_cube))[3]
  
      ;------------------------------[ Begin Thought Process (Bad Pixel Stuff) ]---------------------------------
  
      ; Here we mask out bad pixels
      if do_bad_pixels then begin
         ;print, 'Loading bad pixel mask...'
         ;badpix_mask = '/Users/kevin/Desktop/MWC758/oct_2016_data/copies_data_files/badPixMaskJune2015Data_v1.fits'
         ;badpix_mask = readfits(badpix_mask)
  
         ;Read in our object cube created above from our individual raw FITS files
         obj_cube = readfits(strcompress(output_folder + obj_name + '_cube.fits', /rem))
         restore, strcompress(output_folder + obj_name + '_parang.sav', /rem)
         print, 'Object cube size:', size(obj_cube)
       
         badpix = mean(obj_cube[*,*,0:50], dim=3); mean image for the first 51(?) frames
         badpix_mask = badpix & badpix_mask[*] = 1; Initialize the mask with all 1s
         boxhsize = 8 & threshold = 5 & n_bad_pixels = 0. & n_tests = 0. ; 5 sigma threshold to conclude a pixel is bad
       
         for xx = 1.*boxhsize, (size(obj_cube))[1]-1-(1.*boxhsize) do begin; (size(obj_cube))[2] is our image width (x dimension)
            for yy = 1.*boxhsize,(size(obj_cube))[2]-1-(1.*boxhsize) do begin; (size(obj_cube))[1] is our image height (y dimension)
  	           test = badpix
  	           test[xx,yy] = !values.f_nan
  	           box = test[xx-boxhsize:xx+boxhsize-1, yy-boxhsize:yy+boxhsize-1]
  	           boxstdv = stddev(box, /nan)
  	           this = badpix[xx,yy] - mean(box, /nan)
  	         
  	            ;Do this stuff to find a set a bad pixel
  	           if abs(this) / boxstdv gt threshold then begin
  	              badpix_mask[xx, yy] = 0; 0 means a bad pixel!
  	              ;print, 'Bad pixel found!'
  	              n_bad_pixels += 1.
  	           endif; bad pixel locating if
  	           
  	           n_tests += 1.
  	           
  	           ;Print our progress every 100 frames (every frame is just too jittery) p.s., I should probably make a real widget for this sort of thing :P
  	           if n_tests mod 100 eq 0 then begin
                  cls; clear the IDL terminal
  	              print, 'Testing pixel at (x, y): ', fix(xx), fix(yy)
  	              print, 'n_bad_pixels: ', fix(n_bad_pixels), ' n_tests: ', long(n_tests)
  	              print, 'Progress: ', 100. * n_tests / ((size(obj_cube))[1] * (size(obj_cube))[2]), '%'
  	              print, 'Bad Pixel Percent =', 100. * n_bad_pixels / n_tests, '%'
  	           endif
            endfor; yy for
          
            ;Write the bad pixel mask based on the given y-value that we're on in our loop (why?)
  	        writefits, strcompress(output_folder + obj_name + '_new_badpix_mask.fits', /rem), badpix_mask
         endfor; xx for
  
         ;badpix_mask[where(badpix_mask eq 0.)]=2.
         ;badpix_mask[where(badpix_mask ne 2.)]=0.
         ;badpix_mask[where(badpix_mask eq 2.)]=1.
         ;badpix_mask=fix(badpix_mask)
  
         ; Write the bad pixel mask for a final time since we've looped through every pixel now
         writefits, strcompress(output_folder + obj_name + '_new_badpix_mask.fits', /rem), badpix_mask
  
         n_frames = (size(obj_cube))[3]
         for ii=0,n_frames-1 do begin; go through each frame in the data cube
            if ii mod 10 eq 0 then print, 'Working on frame ', ii, ' / ', n_frames-1
  	        frame = obj_cube[*,*,ii]; get our frame to fix
  	        fixed_frame = frame; initialize a frame that will be our fixed frame soon!
  	        ;frame[where(badpix_mask gt 0)]=!values.f_nan
  	        ;fixpix procedure found in the GPI reduction pipeline files online
  	        fixpix, frame, badpix_mask, fixed_frame, npix=20;, /silent; Fix the frame (why 20?)
  	        obj_cube[*,*,ii] = fixed_frame; reassign our frame in the obj cube to the fixed frame
         endfor; fix frame for
  
         ; write our new FITS with the bad pixels removed
         writefits, strcompress(output_folder + obj_name + '_cube_no_bad_pixels.fits', /rem), obj_cube
  
      endif; bad_pixels if
  
      ;------------------------------[ Begin Thought Process (Sky Subtraction) ]---------------------------------
  
      ; Subtract background stuff from the sky from our newly-created data cube
      if do_sky_subtraction then begin
         obj_cube = readfits(strcompress(output_folder + obj_name + '_cube_no_bad_pixels.fits', /rem))
         ; restores angles and flags from the combining section (why do we need to restore them?) Just if we want to first combine into a cube and then do sky sub later?
         restore, strcompress(output_folder + obj_name + '_parang.sav', /rem)
         ;frames in each nod = 300/bin, e.g. bin=50 => 6 frames in each nod
         frames_per_nod = 200. / coadd
         
         n_frames = (size(obj_cube))[3]
         print, 'Number of frames = ', n_frames
         print, 'Number of nods = ', n_frames / frames_per_nod
  
         firstsky = obj_cube[*,*,16:23] ;manually set for this dataset
         firstsky = mean(firstsky, dim=3)
         writefits,'~/Desktop/firstsky.fits', firstsky
         skys = []
  
  	     sky = firstsky
  
         for ii=0, n_frames-1 do begin
  	        flag_i = flags[ii]
  	        if ii ne n_frames-1 then flag_next_i = flags[ii+1]
  	        print, 'Frame = ', ii, '/', n_frames-1
  	        print, 'This, next flag = ', flag_i, flag_next_i
  	        ;if ii ge 40 then frames_per_nod = 400. / coadd
  	        frame = obj_cube[*,*,ii]
  
  	        skys = [ [[skys]], [[frame]] ]
  	        if flag_next_i ne flag_i then sky = median(skys, dim=3)
  	        if flag_next_i ne flag_i then skys = []
  
            ; This puts the subtraction in sky subtraction!
  	        frame -= sky
  
            ; Assign our frame as the new cube if we're on the first frame, or add it to the array if we're on any subsequent frame
  	        if ii eq 0 then new_cube = frame else new_cube = [ [[new_cube]], [[frame]] ]
         endfor; sky subtraction for
  
         ; write our sky-subtracted "new cube"
         writefits, strcompress(output_folder + obj_name + '_cube_skysub.fits', /rem), new_cube
  
      endif ;sky subtraction if
  
      ;------------------------------[ Begin Thought Process (Split) ]---------------------------------
  
      if do_split then begin
         obj_name = strcompress(obj_name, /rem)
         obj_cube = readfits(output_folder+obj_name + '_cube_skysub.fits')
         restore, filename = output_folder + obj_name + '_parang.sav'; restore angles and flags from the combining step
         oldangles = angles
         ; The third element the size of our cube is the third dimension of the cube, aka the dimension of frames (other two are just x and y for the images)
         n_frames = (size(obj_cube))[3]
         ;frames_per_nod=300./coadd
  
         side = 0
         new_cube = []
         angles = []
  
         ; -1 because our indices start at 0 but our n_frames is just a total number 
         for ii=0, n_frames - 1 do begin
            ;Grab a frame from our cube to work on
  	        newframe = obj_cube[*,*,ii]
        
            ; If we're on any frame other than the first one
            if ii gt 0 then begin
  	           ;if ii mod frames_per_nod eq 0. then print, 'Found nod, switching...'
  	           flag_i = flags[ii]
  	           if ii ne n_frames-1 then flag_next_i = flags[ii+1]; Do this on all but our last frame
  	       
  	           if flag_i ne flag_next_i then begin
  		            if side eq 0 then side = 1 else side=0
  	           endif ; flag_i if
            endif; ii gt 0 if
  	        ;print, ii mod frames_per_nod
  	        ;print, dither_folder, side
  
  	        if dither_folder eq '/dith1/' and side eq 0 then begin
  		         ;print, 'success'
  		         new_cube = [ [[new_cube]], [[newframe]] ]
  		         angles = [angles, oldangles[ii] ]
  	        endif; dith 
  
  	        if dither_folder eq '/dith2/' and side eq 1 then begin
  		         new_cube = [ [[new_cube]], [[newframe]] ]
               angles = [angles, oldangles[ii] ]
  	        endif ;dith2 side = 1 if
         endfor; ii for
  
         if dither_folder eq '/dith1/' then new_cube = new_cube[547-250:547+249, *, *] else new_cube = new_cube[1475-250:1475+249, *, *]
         if output_folder eq left_output_path then new_cube = new_cube[*, 701-250:701+249, *] else new_cube = new_cube[*,259-250:259+249, *] 
  
         ; Write our split stuff, I guess?
         writefits, output_folder + dither_folder + obj_name + '_cube_skysub.fits', new_cube
         save, filename = output_folder + dither_folder + obj_name + '_parang.sav', angles
      endif ;split if
  
      ;------------------------------[ Begin Thought Process (Centering) ]---------------------------------
  
      if do_center then begin
         obj_name = strcompress(obj_name, /rem)
         obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub.fits')
         restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
         obj_cube[where(finite(obj_cube) eq 0)] = 0.
  
         ;correcting bad frames
         ;if output_folder eq '/home/ggdoubleu/OneDrive/Research/HII1348/Testing/processed-left/' then begin
            ;if dither_folder eq '/dith1/' and night eq 2 then obj_cube[*, *, 36] = !values.f_nan
            ;if dither_folder eq '/dith2/' and night eq 2 then obj_cube[*, *, 24:38] = !values.f_nan
         ;endif
  
         ;bad_pixels = [37, 90, 91, 92, 93, 94] - 1
         ;goods = indgen(n_frames)
 
         ;remove, bad_pixels, goods
         ;print, angles
         ;obj_cube = obj_cube[*, *, goods]
         ;angles = angles[goods]
         ;print, angles
  
         ;endif
  
         ;if do_median_subtraction then for ii=0, n_frames-1 do obj_cube[*,*,ii] = obj_cube[*,*,ii] - medframe
  
         ;initialize a new cube to modify with the same data
         obj_cube2 = obj_cube
  
         ; Divide into oblivion if we're blocking a side out
         oblivion = 1000000.
         if do_block_right then obj_cube2[350:499, *, *] /= oblivion
         if do_block_left then obj_cube2[0:150, *, *] /= oblivion
         ;if do_block_left then obj_cube2[700:1023, *, *] /= oblivion
         ;if do_block_right then obj_cube2[0:450, *, *] /= oblivion
         if do_block_bottom then obj_cube2[*, 0:150, *] /= oblivion
         if do_block_top then obj_cube2[*,350:499, *] /= oblivion
         
         ;obj_cube2[*, 0:first_center[0]-25., *] /= oblivion
         ;if do_block_top then obj_cube2[*, first_center[0]+25.:1023., *] /= oblivion
  
         ;obj_cube[*,*,0] = fshift(obj_cube[*,*,0], 300. - first_center[0], 450. - first_center[1])
         ;obj_cube2[*,*,0] = fshift(obj_cube2[*,*,0], 300. -first_center[0], 450. - first_center[1])
  
         obj_cube2[where(finite(obj_cube) eq 0)] = 0.
  
         obj_cube2[where(obj_cube2 lt 0)] /= 100000000. ;reduce negative values so CC doesn't get confused on subtracted PSF.
         half_sizex = 250.
         half_sizey = 250.
  
         n_frames = (size(obj_cube))[3]
         if do_cen_filter then begin
            for ii=0, n_frames - 1 do begin
  	           obj_cube2[*,*,ii] -= smooth(obj_cube2[*,*,ii], 20., /NAN)
  	           obj_cube2[*,*,ii] = smooth(obj_cube2[*,*,ii], 5., /NAN)
            endfor
         endif; do_cen_filter if
  
         ; Returns medframe, whose pixels are the median of those in the 3D obj_cube2
         if do_cen_median_sub then medarr, obj_cube2, medframe
  
         for ii=0, n_frames - 1 do begin
  	        if do_cen_median_sub then obj_cube2[*,*,ii] -= medframe
  	
  	        corr = crosscorr(obj_cube2[*,*,130], obj_cube2[*,*,ii], pmax, /cut)
  	        if pmax(0) ge 0. then dx = half_sizex - pmax(0) else dx = -half_sizex + abs(pmax(0))
  	        if pmax(1) ge 0. then dy = half_sizey - pmax(1) else dy = -half_sizey + abs(pmax(1))
  
  	        print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
  	        obj_cube[*,*,ii] = fshift(obj_cube[*,*,ii], -dx, -dy); fshift shifts frames by non-int values, which we're using to center
         endfor; ii for
  
         ;Basically the same as the medarr used above with obj_cube2
         med = median(obj_cube, dim=3)
         peak = mpfit2dpeak(med, A, /tilt); Gaussian fit
         xx = A[4] & yy = A[5]
     
         ; What is the difference between this fshift and the one above with dx and dy?
         for ii=0, n_frames - 1 do obj_cube[*,*,ii] = fshift(obj_cube[*,*,ii], 250. - xx, 250. - yy)
  
  	     ;obj_cube=obj_cube[*,150:749,*]
  	     ;
         ;Write our centered FITS file
         writefits, output_folder + dither_folder + obj_name + '_cube_skysub_cen.fits', obj_cube
      endif ;center if
  
      ;------------------------------[ Begin Thought Process (Clean) ]---------------------------------
  
      if do_clean then begin
         obj_name = strcompress(obj_name, /rem)
         print, 'Finding bad frames...'
         cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen.fits', hdr)
         restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'

         ;192,195
  
         bad_pixels = []
         goods = findgen((size(cube))[3])
  
         ;initialize a new cube to clean
         cube2 = cube
         ;remove saturated center
         ;for xx = 0, (size(cube))[2] - 1 do for yy = 0, (size(cube))[2] - 1 do if sqrt( (( xx - 250.)^2.) + (( yy - 250.)^2.) ) lt 10 then cube2[xx,yy,*] /= 10000000.
         
         for xx = 0, (size(cube))[2] - 1 do begin
            for yy = 0, (size(cube))[2] - 1 do begin
               ;Pythagorean theorem baby! (everything within r = 20px is divided into oblivion)
               if sqrt( ((xx - 250.)^2.) + ((yy - 250.)^2.) ) gt 20 then cube2[xx,yy,*] /= 10000000.
            endfor
         endfor
         
         ;compute pupil median
         lmed = median(cube2, dim=3)
     
         writefits, '~/Desktop/test-med.fits', lmed
  
         for jj=0, (size(cube))[3] - 1 do begin
  		      maxcor = max(crosscorr(cube2[*,*,jj], lmed))
  		      print, 'Frame ', jj, ' has max correlation value of ', maxcor, ' with pupil median.'
  		      ;First frame vs. the rest
  	        if jj eq 0 then corrs = maxcor else corrs = [corrs, maxcor]
         endfor; max correlation for
  
         bad_pixels = where(corrs lt corr_thresh)
  
         plot, corrs - median(corrs)
         level = corrs
         level[*] = corr_thresh
         ;print, level
         !p.linestyle = 2
         oplot, level - median(corrs)
         !p.linestyle = 0
  
         print, bad_pixels
         print, 'Found ', n_elements(bad_pixels), ' bad frames.'
         print, 'Average correlation = ', mean(corrs)
         print, 'Median correlation = ', median(corrs)
         print, 'Standard deviation = ', stdev(corrs)
         ;hak
         print, 'Bad frame percentage = ', 100. * n_elements(bad_pixels) / n_elements(corrs), '%'
   
         if n_elements(bad_pixels) ge 1 and bad_pixels ne -1 then begin
           remove, bad_pixels, goods, angles
           left_bad_pixels_cube = cube[*,*,bad_pixels]
           writefits, '~/Desktop/test-bad_pixels.fits', left_bad_pixels_cube
         endif
  
         cube = cube[*,*,goods]; Reassign our cube to only have our good pixels
  	     medarr, cube, medframe; Get a median frame for our new bad-pixel-free cube
  
         ;Write our clean cube, our bad pixels, and our median frame
  	     writefits, output_folder + dither_folder + obj_name + '_cube_skysub_cen_clean.fits', cube, hdr
  	     writefits, output_folder + dither_folder + obj_name + '_cube_skysub_cen_bad_pixels.fits', left_bad_pixels_cube, hdr
  	     writefits, output_folder + dither_folder + obj_name + '_pupil.fits', medframe
  	     save, filename = output_folder + dither_folder + obj_name + '_parang_clean.sav', angles
  
         cgcleanup; Found online, apparently it closes any graphics or widget windows? Did we create any??
      endif; clean if
  
      ;------------------------------[ Begin Thought Process (Planet Injection) ]---------------------------------
  
      if do_inject then begin
         obj_name = strcompress(obj_name, /rem)
         obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen_clean.fits')
         ref_file = output_folder + dither_folder + obj_name + '_pupil.fits' ;unsaturated data can just use the pupil image
         ref = readfits(ref_file)
  
         restore, filename = output_folder + dither_folder + obj_name + '_parang_clean.sav'
  	     
  	     if runs eq 3 or runs eq 4 then truenorth = 0.59 else truenorth = -1.39
  
         ;derotate
         n_frames = (size(obj_cube))[3]
         for ii=0, n_frames - 1 do begin
  	        print, 'Derotating by ', -angles[ii] - truenorth
  	        obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii] - truenorth, /interp)
         endfor; derotate for
  
         big_ref = obj_cube[*,*,0]
         big_ref[*] = 0
  
         ;build bigger reference image
         for xx=0, (size(ref))[1] - 1 do begin
  	        for yy=0, (size(ref))[1] - 1 do begin
  		         big_ref[250. - (size(ref))[1] /2. + xx, 250. - (size(ref))[1] / 2. + yy] = ref[xx,yy]
  	        endfor;yy for
         endfor;; xx for
  
         ;inject planets
         for ii=0, (size(planet_r))[1] - 1 do begin
            print, 'Injecting planet', ii
  	        big_ref_ii = big_ref * planet_contrast[ii]
  	
  	        xshift = planet_r[ii] * (1. / pixel_scale) * Cos(planet_theta[ii])
  	        yshift = planet_r[ii] * (1. / pixel_scale) * Sin(planet_theta[ii])
  	        big_ref_ii = fshift(big_ref_ii, xshift, yshift)
  	
            for jj=0, n_frames - 1 do obj_cube[*,*,jj] += big_ref_ii
         endfor; planet injection for
  
         ;rotate back to pupil stabilized orientation
         for ii=0, n_frames - 1 do begin
  	        print, 'Rotating by ', angles[ii]
  	        obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], angles[ii] + truenorth, /interp)
         endfor;rotate back for
  
  	     writefits, output_folder + dither_folder + obj_name + '_cube_skysub_cen_inj.fits', obj_cube
      endif ;inject if
  
      ;------------------------------[ Begin Thought Process (Rotate) ]---------------------------------
  
      if do_rotate then begin
         obj_name = strcompress(obj_name, /rem)
         obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen_clean.fits')
         if use_injection then obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen_inj.fits')
  
         n_frames = (size(obj_cube))[3] 
         if do_destripe then begin
            print, 'destriping 90 degrees...'
            for ii=0, n_frames - 1 do obj_cube[*,*,ii] = destripe(obj_cube[*,*,ii], 90., clip_level=0.0, /nodisp)
            print, 'destriping 0 degrees...'
            for ii=0, n_frames - 1 do obj_cube[*,*,ii] = destripe(obj_cube[*,*,ii], 0., clip_level=0.0, /nodisp)
         endif; destripe if
  
         if filter gt 1 then for iii=0, n_frames - 1 do obj_cube[*,*,iii] -= smooth(obj_cube[*,*,iii], filter)
      
         restore, filename = output_folder + dither_folder + obj_name + '_parang_clean.sav'
      
      	 medarr, obj_cube, medframe; medframe is the output here of calling medarr on our obj_cube
      	 
      	 ;Initialize our cube for Angular Differential Imaging
      	 adi_cube = obj_cube
      	
      	 if runs eq 3 or runs eq 4 then truenorth = 0.59 else truenorth = -1.39
      
         for ii=0, n_frames - 1 do begin
      	    print, 'Rotating by ', angles[ii]
      	    obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii] - truenorth, /interp)
      	    adi_cube[*,*,ii] -= medframe
      	    adi_cube[*,*,ii] = rot(adi_cube[*,*,ii], -angles[ii] - truenorth, /interp)
         endfor; rotate if
      
      	 writefits, output_folder + dither_folder + obj_name + '_cube_skysub_cen_filt_derot.fits', obj_cube
      	
      	 if combine_type eq 'median' then medarr, obj_cube, medframe
      	 if combine_type eq 'mean' then medframe = mean(obj_cube, dim=3)
      	 if combine_type eq 'nwadi' then medframe = nw_ang_comb(obj_cube, angles)
      	
      	 writefits, strcompress(output_folder + dither_folder + obj_name + '_median_derot.fits', /rem), medframe
      
      	 if combine_type eq 'median' then medarr, adi_cube, adiframe
      	 if combine_type eq 'mean' then adiframe = mean(adi_cube, dim=3)
      	 if combine_type eq 'nwadi' then adiframe = nw_ang_comb(adi_cube, angles)
      
      	 adiframe[where(finite(adiframe) ne 1)] = 0.
      	
      	 writefits, strcompress(output_folder + dither_folder + obj_name + '_median_derot_adi.fits', /rem), adiframe
      	
         size = 500.
      	 width = 8.72059;(3.8*1E-6) / (8.4) * 206265. / 0.0107; Where does this come from? 
      	 print, 'PSF Width: ', width
      	 PSF = psf_Gaussian(npixel=size, FWHM=[width, width])
      	 PSFN = PSF / MAX(PSF); N for normalized?
      	 adiframe_c = convolve(adiframe, PSFN)
      	
      	 writefits, strcompress(output_folder + dither_folder + obj_name + '_median_derot_adi_conv.fits', /rem), adiframe_c
      
      endif ;rotate if
      
      ;------------------------------[ Begin Thought Process (KLIP) ]---------------------------------
      
      if do_klip then begin
         obj_name = strcompress(obj_name, /rem)
         obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen_clean.fits')
         if use_injection then obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen_inj.fits')
      
         n_frames = (size(obj_cube))[3]
         if do_destripe then begin
            print, 'destriping 90 degrees...'
            for ii=0, n_frames - 1 do obj_cube[*,*,ii] = destripe(obj_cube[*,*,ii], 90., clip_level=0.0, /nodisp)
            print, 'destriping 0 degrees...'
            for ii=0, n_frames - 1 do obj_cube[*,*,ii] = destripe(obj_cube[*,*,ii], 0., clip_level = 0.0, /nodisp)
         endif; destripe if
      
         if filter gt 1 then for iii=0, n_frames - 1 do obj_cube[*,*,iii] -= smooth(obj_cube[*,*,iii], filter)
      
         restore, filename = output_folder + dither_folder + obj_name + '_parang_clean.sav'
      
         if bin gt 1 and bin_type eq 'median' then begin
      	    st=1
      	    
      	    for ii=0., n_frames - 1 do begin
      		     print, fix(ii + 1.) mod fix(bin)
      		     
      		     if fix(ii + 1.) mod fix(bin) eq 1 then binned = obj_cube[*,*,ii] else binned = [ [[binned]], [[ obj_cube[*,*,ii] ]] ]
      		     
      		     if fix(ii + 1.) mod fix(bin) eq 0 then begin
      			      print,'Binning left frames...'
      			      
      			      if st then begin 
      			         ;medarr,binned,binned
      				       binned = median(binned, dim=3, /even)
      			         binned_cube = binned 
      			         st = 0 
      			         
      			      endif else begin
      			         ;medarr,binned,binned
      				       binned = median(binned, dim=3, /even)
      			         binned_cube = [[[binned_cube]], [[binned]]]
      			      endelse
      			      
      			      print, size(binned_cube)
      		     endif; fix mod eq 0 if
      	    endfor; ii for
      	    
      	    print, size(binned_cube)
      	    st=1
      	    
      	    for ii=0, n_frames - 1. do begin
      		     if fix(ii + 1.) mod fix(bin) eq 1 then binned_angle = angles[ii] else binned_angle += angles[ii]
      		     
      		     if fix(ii + 1.) mod fix(bin) eq 0 then begin
      			      if st then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles = [[[binned_angles]], [[binned_angle / bin]]]
      		     endif
      	    endfor; ii for	
      	
      	    obj_cube = binned_cube
      	    angles = binned_angles
      	       
         endif; bin gt 1 if
      
         if bin gt 1 and bin_type eq 'mean' then begin
      	    st=1
      	    
      	    for ii=0., n_frames - 1. do begin
      		     print, fix(ii + 1.) mod fix(bin)
      		     if fix(ii + 1.) mod fix(bin) eq 1 then binned = obj_cube[*,*,ii] else binned += obj_cube[*,*,ii] 
      		     if fix(ii + 1.) mod fix(bin) eq 0 then begin
      			      print,'Binning left frames...'
      			
      			      if st then begin 
      			         ;medarr,binned,binned
      			         binned_cube = binned / float(bin)
      			         st=0 
      			      endif else begin
      			         ; medarr,binned,binned
      			         binned_cube = [ [[binned_cube]], [[binned / float(bin)]] ]
      			      endelse
      			      
      			      print, size(binned_cube)
      			      
      		     endif; fix mod eq 0 if
      	    endfor; ii for
      	    
      	    print, size(binned_cube)
      	    st=1
      	    
      	    for ii=0, n_frames - 1. do begin
      		     if fix(ii + 1.) mod fix(bin) eq 1 then binned_angle = angles[ii] else binned_angle += angles[ii]
      		     
      		     if fix(ii + 1.) mod fix(bin) eq 0 then begin
      			      if st then begin binned_angles = binned_angle / bin & st=0 & endif else binned_angles = [ [[binned_angles]], [[binned_angle / bin]] ]
      		     endif
      	    endfor; ii for
      	
      	    obj_cube = binned_cube
      	    angles = binned_angles
      	    
         endif; bin gt 1 if
      
         lbc = 0
         bad_pixels = fltarr(n_frames)
         bad_pixels[*] = 0.
      
         obj_cube[where(finite(obj_cube) eq 0)] = 0.
      
         for ii=0, n_frames - 1 do begin
      	    if total(obj_cube[*,*,ii]) eq 0 then begin
      	       print, ii
      		     remove, ii - lbc, angles
      		     bad_pixels[ii] = 1.
      		     lbc += 1
      	    endif; total eq 0 if
      	    ;hsize = h_size / 2.
      	    ;center with centroid
      	    ;cntrd, left[*,*,ii], hsize - 1. , hsize - 1. , xcc, ycc, 5
      	    ;left[*,*,ii] = fshift(left[*,*,ii], -(xcc-(hsize-1.)), -(ycc-(hsize-1.)))	
         endfor; ii for
         
         obj_cube = obj_cube[*, *, where(bad_pixels eq 0.)]
      
         if klip_fraction then begin
      	    obj_cube = obj_cube[*,*, (start_frame / 100.) * (n_frames):(end_frame / 100.) * (n_frames - 1.)]
      	    angles = angles[(start_frame / 100.) * (n_frames):(end_frame / 100.) * (n_frames - 1.)]
         endif
      
         klip_cube = obj_cube
      
         num = 1 ;vestigal
      
         ; Where do I get the adiklip function for this?
         n_frames = (size(obj_cube))[3]
         for ii=0, n_frames - 1 do begin
            print, '-------[ KLIPing image ', ii, ' out of ', n_frames - 1, ' ]----------'
            if not do_hyper and do_annmode ne 1 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=angles, wl=4.7,diam=8.4, pixelscale=0.0107, ang_sep=ang_sep,max_angle=max_angle, obj_name=obj_name,n_rings=n_rings, wr =wr, n_angle =n_angle, num=758) 
            if do_hyper and do_annmode ne 1 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=angles, wl=4.7,diam=8.4, pixelscale=0.0107, ang_sep=ang_sep,max_angle=max_angle, obj_name=obj_name,n_rings=n_rings, wr =wr, n_angle =n_angle, num=758, /hyper) 
            if do_annmode then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=angles, wl=4.7,diam=8.4, pixelscale=0.0107, ang_sep=ang_sep, max_angle=max_angle, obj_name=obj_name,n_rings=n_rings, wr =wr, n_angle =n_angle, num=758, annmode_inout=annmode_inout) 
         endfor; klip ii for
      
      	 truenorth = -1.39
      	 if runs eq 3 or runs eq 4 then truenorth = 0.59
      
         for ii=0, n_frames - 1 do begin
      	    print, 'Rotating by ', angles[ii]
      	    klip_cube[*,*,ii] = rot(klip_cube[*,*,ii], -angles[ii] - truenorth, /interp)
      	    obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii] - truenorth, /interp)
      
      	    framei = klip_cube[*,*,ii]
      	    frameifull = obj_cube[*,*,ii]
      	
      	    ;fill in the rest of the image
      	    if fill then framei[where(finite(framei) eq 0)] = frameifull[where(finite(framei) eq 0)]
      	    klip_cube[*,*,ii] = framei
         endfor; rotating ii for
      
      	 writefits, output_folder + dither_folder + obj_name + '_cube_klip.fits', klip_cube
      	
      	 if combine_type eq 'median' then medarr, klip_cube, medframe
      	 if combine_type eq 'mean' then medframe = mean(klip_cube, dim=3)
      	 if combine_type eq 'nwadi' then medframe = nw_ang_combine(klip_cube, angles)
      	
      	 writefits, strcompress(output_folder + dither_folder + '/klip/' + obj_name + '_median_klip' + suffix + '.fits',/rem), medframe
      	
      	 size = 500.
      	 width = 8.72059 ;(3.8*1E-6) / (8.4) * 206265. / 0.0107; What are these numbers?
      	 print, 'PSF Width: ', width
      	 PSF = psf_Gaussian(npixel=size, FWHM=[width,width])
      	 PSFN = PSF / MAX(PSF)
      
      	 medframe[where(finite(medframe) ne 1)] = 0.
      	 medframe_c = convolve(medframe, PSFN)
      			
      	 writefits, strcompress(output_folder + dither_folder + '/klip/' + obj_name + '_median_klip_conv' + suffix + '.fits', /rem), medframe_c
      
         klipframe = medframe
      endif ;klip if
      
      ;------------------------------[ Begin Thought Process (RDI KLIP) ]---------------------------------
      
      if do_rdi_klip then begin
         obj_cube = readfits(output_folder + dither_folder + obj_name + '_cube_skysub_cen_filt.fits')
         ref_cube = readfits(ref_folder + dither_folder + ref_name + '_cube_skysub_cen_filt.fits')
         
         n_frames = (size(obj_cube))[3]
         if filter gt 1 then for iii=0, n_frames - 1 do obj_cube[*,*,iii] -= smooth(obj_cube[*,*,iii], filter)
      
         restore, filename = output_folder + dither_folder + obj_name + '_parang.sav'
      
         if klip_fraction then begin
      	    obj_cube = obj_cube[*,*, start_frame:end_frame]
      	    angles = angles[start_frame:end_frame]
         endif; klip fraction if
      	
         if rdi_high_pass gt 0 then begin
            for ii=0, n_frames - 1 do begin
      		     obj_cube[*,*,ii] -= smooth(obj_cube[*,*,ii], rdi_high_pass, /NAN)
      			endfor
      			
      			for ii=0, (size(ref_cube))[3] - 1 do begin
      			   ref_cube[*,*,ii] -= smooth(ref_cube[*,*,ii],rdi_high_pass,/NAN)
      			endfor
         endif; RDI high pass if
      
      	 if do_cross_corr_sci then begin
      	    hs = 100.
      	    subx = 512.
      	    suby = 512.
      	    subcube_obj = obj_cube[subx-hs:subx+hs-1., suby-hs:suby+hs-1., *]
      	    subcube_ref = ref_cube[subx-hs:subx+hs-1., suby-hs:suby+hs-1., *]
      
      	    half_size = hs
      					
      	    for ii=0, n_frames - 1 do begin
      		     subcube_obj[*,*,ii] -= smooth(subcube_obj[*,*,ii], 15.)
      		     subcube_obj[*,*,ii] = smooth(subcube_obj[*,*,ii], 5.)
      	    endfor
      
            for ii=0, (size(ref_cube))[3] - 1 do begin
      		     subcube_ref[*,*,ii] -= smooth(subcube_ref[*,*,ii], 15.)
      	       subcube_ref[*,*,ii] = smooth(subcube_ref[*,*,ii], 5.)
      	
      		     corr = crosscorr(subcube_obj[*,*,0], subcube_ref[*,*,ii], pmax, /cut)
      		     if pmax(0) ge 0. then dx = half_size - pmax(0) else dx = -half_size + abs(pmax(0))
      		     if pmax(1) ge 0. then dy = half_size - pmax(1) else dy = -half_size + abs(pmax(1))
      
      	       print, 'Shifting frame', ii, ' / ', n_frames - 1, ' by: ', String(-dx), String(-dy)
      		     ref_cube[*,*,ii] = fshift(ref_cube[*,*,ii], -dx, -dy)
            endfor
      
      		  medarr, obj_cube, medframe
      		  medarr, ref_cube, medframeref
      
      		  writefits,output_folder+dither_folder+obj_name+'_rdi_align_med.fits', [ [[medframe]], [[medframeref]] ]
      	 endif; do_cross_corr_sci if
      
         ;initialize a new cube for KLIP
         klip_cube = obj_cube
         num = 1 ;vestigal
      
         for ii = 0,n_frames - 1 do begin
            print, '-------[ KLIPing image ', ii, ' out of ', n_frames - 1, ' ]----------'
            if not do_hyper then klip_cube[*,*,ii] = rdiklip(obj_cube,ref_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=angles, wl=3.8,diam=8.4, pixelscale=0.0107, ang_sep=ang_sep, obj_name=obj_name,n_rings=n_rings, wr =wr, n_angle =n_angle) 
            if do_hyper then klip_cube[*,*,ii] = rdiklip(obj_cube,ref_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, posang=angles, wl=3.8,diam=8.4, pixelscale=0.0107, ang_sep=ang_sep, obj_name=obj_name,n_rings=n_rings, wr =wr, n_angle =n_angle, /hyper)
         endfor; KLIP for
      
      	 if runs eq 3 or runs eq 4 then truenorth = 0.59 else truenorth = -1.39
      	 
         for ii = 0, n_frames - 1 do begin
            print, 'Rotating by ', angles[ii]
            klip_cube[*,*,ii] = rot(klip_cube[*,*,ii], -angles[ii] - truenorth, /interp)
            obj_cube[*,*,ii] = rot(obj_cube[*,*,ii], -angles[ii] - truenorth, /interp)
      
            framei = klip_cube[*,*,ii]
            frameifull = obj_cube[*,*,ii]
      	
            ;fill in the rest of the image
            if fill then framei[where(finite(framei) eq 0)] = frameifull[where(finite(framei) eq 0)]
            klip_cube[*,*,ii] = framei
         endfor; Rotate for
      
         writefits, output_folder + dither_folder + obj_name + '_cube_rdiklip.fits', klip_cube
         medarr, klip_cube, medframe
         writefits, output_folder + dither_folder + '/klip/' + obj_name + '_median_rdiklip' + suffix + '.fits', medframe
      endif ;do rdiklip if
  
      ;------------------------------[ Begin Thought Process (Finishing Touches per run) ]---------------------------------
  
      if do_klip then begin
         if runs eq 1 then nods = klipframe else nods = [[[nods]], [[klipframe]]]
         writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_nods_klip' + suffix + '.fits', /rem), nods
      endif
  
      if do_rotate then begin
         if runs eq 1 then adinods = adiframe else adinods = [[[adinods]], [[adiframe]]]
         writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_nods_adi' + suffix + '.fits', /rem), adinods
      endif
  
      ;The final run! Let's write a bunch of FITS!
      if runs eq 4 then begin
         if do_klip then begin
            e = mean(nods, dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_total_klip' + suffix + '.fits', /rem), e
            
            e = mean(nods[*,*,0:1], dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_left_klip' + suffix + '.fits', /rem), e
            
            e = mean(nods[*,*,2:3], dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_right_klip' + suffix + '.fits', /rem), e
            
            ;used to combine injections from both nights for improved SNR. Manually set to inj for now.
            e = mean(nods, dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_total_klip.fits', /rem), e
            
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_klip_nod1.fits', /rem), nods[*,*,0]
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_klip_nod2.fits', /rem), nods[*,*,1]
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_klip_nod3.fits', /rem), nods[*,*,2]
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_klip_nod4.fits', /rem), nods[*,*,3]
            ;Where does this correction factor come from?
            find_sources, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_total_klip.fits', /rem), reference=output_folder+dither_folder+obj_name+'_pupil.fits', platescale=0.0107, correction_factor=2.8628373, fwhm=8.7
            
            size = 500.
            width = 10.7860;(4.7*1E-6) / (8.4) * 206265. / 0.0107
            print, 'PSF Width: ',width
            PSF = psf_Gaussian(npixel=size, FWHM=[width, width])
            PSFN = PSF / MAX(PSF)
            ec = convolve(e, PSFN)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_total_klip_conv.fits', /rem), ec
         endif; do_klip if
         
         ;if runs eq 2 or runs eq 4 then if do_klip then	e=median(nods, dim=3, /even)	
         ;if runs eq 3 then if do_klip then	e=median(nods, dim=3)	
         ;if do_klip then	writefits,strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/'+obj_name+'_klip_mednods'+suffix+'.fits',/rem), e
  
         if do_rotate then begin
            e = mean(adinods, dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_total_adi.fits', /rem), e
            ;What is the correction factor for?
            find_sources, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_total_adi.fits', /rem), reference= ref_file, platescale=0.0107, correction_factor=3.36813, outrad=150, inrad=annmode_inout[0], fwhm=8.7
         
            e = mean(adinods[*,*,0:1], dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_left_adi' + suffix + '.fits', /rem), e
            
            e = mean(adinods[*,*,2:3], dim=3)
            writefits, strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/' + obj_name + '_right_adi' + suffix + '.fits', /rem), e
         endif; do_rotate if
  
         ;if runs eq 2 or runs eq 4 then if do_rotate then 	e=median(adinods, dim=3, /even)	
         ;if  runs eq 3 then if do_rotate then 	e=median(adinods, dim=3)	
         ;if do_rotate then 	writefits,strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/'+obj_name+'_adi_mednods'+suffix+'.fits',/rem), e
         ;if do_klip then	find_sources,strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/'+obj_name+'_klip_nod1.fits',/rem),reference=output_folder+dither_folder+obj_name+'_pupil.fits',platescale=0.0107,correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),outrad=annmode_inout[1]-20.,inrad=annmode_inout[0],fwhm=8.7
         ;if do_klip then	find_sources,strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/'+obj_name+'_klip_nod2.fits',/rem),reference=output_folder+dither_folder+obj_name+'_pupil.fits',platescale=0.0107,correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),outrad=annmode_inout[1]-20.,inrad=annmode_inout[0],fwhm=8.7
         ;if do_klip then	find_sources,strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/'+obj_name+'_klip_nod3.fits',/rem),reference=output_folder+dither_folder+obj_name+'_pupil.fits',platescale=0.0107,correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),outrad=annmode_inout[1]-20.,inrad=annmode_inout[0],fwhm=8.7
         ;if do_klip then	find_sources,strcompress('/home/ggdoubleu/OneDrive/Research/HII1348/Testing/combined/'+obj_name+'_klip_nod4.fits',/rem),reference=output_folder+dither_folder+obj_name+'_pupil.fits',platescale=0.0107,correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),outrad=annmode_inout[1]-20.,inrad=annmode_inout[0],fwhm=8.7
         ;if do_klip then f=readfits('/Volumes/Storage/LBT/MWC758/combined/MWC758_total_klip.fits')
         ;if do_klip then begin
  	     ;for xx=0., 599. do for yy=0., 599. do if sqrt( (xx-300.)^2. + (yy-300.)^2. ) lt 25 then e[xx,yy]=0.
  	     ;for xx=0., 599. do for yy=0., 599. do if sqrt( (xx-300.)^2. + (yy-300.)^2. ) lt 25 then f[xx,yy]=0.	
         ;endif
         ;writefits,'~/Desktop/LBTI_MWC758_Lp_combined.fits',[[[(e+f)/2.]],[[f]],[[e]]]
      endif; runs eq 4 if
      
   endfor; runs = 1,4 for
   
   ;------------------------------[ Begin Thought Process (Actual Finishing Touches) ]---------------------------------
   
   print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 1440., ' minutes.'; 1440 min / day (Julian dates are measured in days)
  
   ;combine_mwc758
   ;endfor
   ;endfor
   ;endfor
   ;endfor
   ;endfor
   ;endfor
  
   ; Why the two identical prints?
   ;print, 'Completed reduction in ', (systime(/JULIAN) - start_time) * 86400. / 60.,' minutes.'
end; That's all, folks!
