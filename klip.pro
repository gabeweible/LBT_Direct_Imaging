pro klip, obj, cube_folder, use_injection, do_destripe, filter, bin, bin_type,$
 do_hyper, do_annmode, combine_type, klip_fraction, start_frame, end_frame, fill,$
  k_klip, angsep, anglemax, nrings, wr, n_ang, annmode_inout, suffix, ct,$
   do_cen_filter, coadd, rho=rho, theta=theta, contrast=contrast, trial=trial,$
    fs=fs, neg_inj=neg_inj
newline = string(10B)

for runs=1,4 do begin

   ; Do this for runs eq 1 and runs eq 3
   if runs mod 2 then dither_folder = 'dith1/' else dither_folder = 'dith2/'
   ;Do this for runs eq 1 and runs eq 2
   if runs lt 3 then begin
      output_folder = cube_folder + 'processed_left/'
      truenorth = -1.39
   endif else begin
      truenorth = 0.59
      output_folder = cube_folder + 'processed_right/'
   endelse

obj=strcompress(obj,/rem)

if use_injection then obj_cube = readfits(output_folder + dither_folder + obj + string(ct) +  '_cube_skysub_cen_clean_inj.fits') else obj_cube = readfits(output_folder + dither_folder + obj + string(ct) +  '_cube_skysub_cen_clean.fits')

if do_destripe then begin
   print, 'destriping 90 degrees...'
   for ii=0, (size(obj_cube))(3)-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],90.,clip_level=0.0,/nodisp)
   print, 'destriping 0 degrees...'
   for ii=0, (size(obj_cube))(3)-1 do obj_cube[*,*,ii]=destripe(obj_cube[*,*,ii],0.,clip_level=0.0,/nodisp)
endif


if filter gt 1 then for iii=0,(size(obj_cube))(3)-1 do obj_cube[*,*,iii]=obj_cube[*,*,iii]-smooth(obj_cube[*,*,iii],filter)

restore,filename=output_folder + dither_folder + obj + string(ct) +  '_parang_clean.sav'

if bin gt 1 and bin_type eq 'median' then begin
   st=1
   for ii=0.,(size(obj_cube))(3)-1. do begin
      print, fix(ii+1.) mod fix(bin)
      if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=[ [[binned]],[[ obj_cube[*,*,ii] ]] ]
      if fix(ii+1.) mod fix(bin) eq 0 then begin
      
         print,'Binning left frames...'
         
         if st eq 1 then begin 
         
         ;medarr,binned,binned

            binned=median(binned,dim=3,/even)
         binned_cube = binned 
         st=0 
         
          endif else begin
          ;medarr,binned,binned

            binned=median(binned,dim=3,/even)

          binned_cube=[[[binned_cube]],[[binned]]]
          
          endelse
         print, size(binned_cube)
      endif
   endfor   
   print, size(binned_cube)
   st=1
   for ii=0,(size(obj_cube))(3)-1. do begin
      if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
      if fix(ii+1.) mod fix(bin) eq 0 then begin
         if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin]]]
      endif
   endfor   
   
   obj_cube=binned_cube
   angles=binned_angles
endif

if bin gt 1 and bin_type eq 'mean' then begin
   st=1
   for ii=0.,(size(obj_cube))(3)-1. do begin
      print, fix(ii+1.) mod fix(bin)
      if fix(ii+1.) mod fix(bin) eq 1 then binned=obj_cube[*,*,ii] else binned=binned + obj_cube[*,*,ii] 
      if fix(ii+1.) mod fix(bin) eq 0 then begin
      
         print,'Binning left frames...'
         
         if st eq 1 then begin 
         
         ;medarr,binned,binned
         binned_cube = binned / float(bin)
         st=0 
         
          endif else begin
         ; medarr,binned,binned

          binned_cube=[[[binned_cube]],[[binned/float(bin)]]]
          
          endelse
         print, size(binned_cube)
      endif
   endfor   
   print, size(binned_cube)
   st=1
   for ii=0,(size(obj_cube))(3)-1. do begin
      if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=binned_angle+angles[ii]
      if fix(ii+1.) mod fix(bin) eq 0 then begin
         if st eq 1 then begin binned_angles = binned_angle/bin & st=0 & endif else binned_angles=[[[binned_angles]],[[binned_angle/bin]]]
      endif
   endfor   
   
   obj_cube=binned_cube
   angles=binned_angles
endif


lbc=0
bads=fltarr( (size(obj_cube))(3) )
bads[*]=0.

obj_cube[where( finite( obj_cube ) eq 0 )]=0.

for ii=0, (size(obj_cube))(3)-1 do begin
   
   if total(obj_cube[*,*,ii]) eq 0 then begin
      print, ii
      remove, ii-lbc, angles
      bads[ii]=1.
      lbc=lbc+1
   endif
   
   ;hsize=h_size/2.
   ;center with centroid
   ;cntrd, left[*,*,ii],hsize-1.,hsize-1.,xcc,ycc, 5
   ;left[*,*,ii]=fshift(left[*,*,ii],-(xcc-(hsize-1.)),-(ycc-(hsize-1.)))  
endfor
obj_cube=obj_cube[*,*,where(bads eq 0.)]



if klip_fraction eq 1 then begin
   obj_cube=obj_cube[*,*,(start_frame/100.)*((size(obj_cube))(3)):(end_frame/100.)*((size(obj_cube))(3)-1.)]
   angles=angles[(start_frame/100.)*((size(obj_cube))(3)):(end_frame/100.)*((size(obj_cube))(3)-1.)]
endif


klip_cube=obj_cube


num=1 ;vestigal

for ii=0, (size(obj_cube))(3)-1 do begin

Print, '-------[ KLIPing image ', ii, ' out of ', (size(obj_cube))(3)-1, ' ]----------'
print, 'Run:', runs
if keyword_set(rho) then begin
   print, newline, 'Separation:', rho
   print, 'Theta:', theta
   print, 'Contrast:', contrast, newline
endif

if do_hyper ne 1 and do_annmode ne 1 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=4.7,diam=8.4, pixelscale=0.0107, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, num=758) 

if do_hyper eq 1 and do_annmode ne 1 then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, trial=trial, anglemask=anglemask, distmask=distmask, posang=angles, wl=4.7,diam=8.4, pixelscale=0.0107, angsep=angsep,anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, num=758, /hyper) 

if do_annmode then klip_cube[*,*,ii] = adiklip(obj_cube, k_klip, target=ii, anglemask=anglemask, distmask=distmask, trial=trial, posang=angles, wl=4.7,diam=8.4, pixelscale=0.0107, angsep=angsep, anglemax=anglemax, obj=obj,nrings=nrings, wr =wr, n_ang =n_ang, num=758, annmode_inout=annmode_inout) 


endfor

   truenorth=-1.39
   if runs eq 3 or runs eq 4 then truenorth=0.59

for ii=0, (size(obj_cube))(3)-1 do begin
   print, 'Rotating by ', angles[ii]
   klip_cube[*,*,ii]=rot(klip_cube[*,*,ii],-angles[ii]-truenorth,/interp)
   obj_cube[*,*,ii]=rot(obj_cube[*,*,ii],-angles[ii]-truenorth,/interp)

   framei=klip_cube[*,*,ii]
   frameifull=obj_cube[*,*,ii]
   
   ;fill in the rest of the image
   if fill eq 1 then framei[where(finite(framei) eq 0)]=frameifull[where(finite(framei) eq 0)]
   klip_cube[*,*,ii]=framei
endfor

   writefits,output_folder+dither_folder+obj+ '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '_cube_klip.fits', klip_cube
   
   if combine_type eq 'median' then medarr, klip_cube, medframe
   
   if combine_type eq 'mean' then medframe=mean(klip_cube,dim=3)
   if combine_type eq 'nwadi' then medframe=nw_ang_comb(klip_cube,angles)
   
   writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip'+suffix+ '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe
   
      sz=500.
         width=(3.8*1E-6) / (8.4) * 206265. / 0.0107
         print, 'PSF Width: ',width
         PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
         PSFN = PSF/MAX(PSF)

         medframe[where(finite(medframe) ne 1)]=0.
         medframe_c = convolve(medframe, PSFN)
         
         writefits,strcompress(output_folder+dither_folder+'/klip/'+obj+'_median_klip_conv'+suffix+ '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits',/rem), medframe_c

klipframe=medframe
   
   if runs eq 1 then nods = klipframe else nods = [[[nods]], [[klipframe]]]
   writefits, strcompress(cube_folder+ 'combined/' + obj + '_nods_klip' + suffix + '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +  '.fits', /rem), nods
   
   if runs eq 4 then begin
      e = mean(nods, dim=3)
      writefits, strcompress(cube_folder + 'combined/' + obj + '_total_klip' + suffix + '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), e

      e = mean(nods[*,*,0:1], dim=3)
      writefits, strcompress(cube_folder + 'combined/' + obj + '_left_klip' + suffix + '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), e

      e = mean(nods[*,*,2:3], dim=3)
      writefits, strcompress(cube_folder + 'combined/' + obj + '_right_klip' + suffix + '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem), e

      ;used to combine injections from both nights for improved SNR. Manually set to inj for now.
      e = mean(nods, dim=3)

      super_suffix = cube_folder + 'combined/' + obj + '_bin_' + string(sigfig(bin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+string(sigfig(k_klip,1))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj)
      if keyword_set(trial) then super_suffix += '_trial_' + string(sigfig(trial, 4))

      writefits, strcompress(super_suffix +   '_total_klip.fits', /rem), e

      writefits, strcompress(super_suffix +   '_klip_nod1.fits', /rem), nods[*,*,0]
      writefits, strcompress(super_suffix +   '_klip_nod2.fits', /rem), nods[*,*,1]
      writefits, strcompress(super_suffix +   '_klip_nod3.fits', /rem), nods[*,*,2]
      writefits, strcompress(super_suffix +   '_klip_nod4.fits', /rem), nods[*,*,3]
      ;Where does this correction factor come from?
      ; I'm having trouble with trying to *not* manually type in the folder here.
      if keyword_set(fs) then begin
	print, 'FS = ', fs
      	if fs eq 1 then begin
         find_sources,strcompress(super_suffix +   '_total_klip.fits',/rem),reference=output_folder+dither_folder+obj+string(ct) + '_pupil.fits',platescale=0.0107,correction_factor=((2.5E-4)/0.00013041987)*((2.0E-4)/0.00013391511),fwhm=8.7,ct=ct,filter=filter,do_cen_filter=do_cen_filter
        endif; find_sources if
      endif; keyword_set(fs) if
      size = 500.
      width = 10.7860;(4.7*1E-6) / (8.4) * 206265. / 0.0107
      print, 'PSF Width: ',width
      PSF = psf_Gaussian(npixel=size, FWHM=[width, width])
      PSFN = PSF / MAX(PSF)
      ec = convolve(e, PSFN)
      writefits, strcompress(super_suffix +  '_total_klip_conv.fits', /rem), ec
   endif
 
endfor; Runs for   
   
end; end procedure
