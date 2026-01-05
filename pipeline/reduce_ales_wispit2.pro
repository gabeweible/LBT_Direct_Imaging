pro reduce_ales_wispit2


;old=readfits('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/reduced/PDS201_sdiklip.fits')

folder='/Volumes/RAID36TB/LMIRCam/TYC5709_2/reduced/'
obj='TYC5709'
pixelscale=0.0345


;implement preprocessing here



do_make_skys=0

do_sub_sky=0
	no_sky=0

do_preproc=0 ;run Jordan's pre-processing script, does CDS here now

do_sort=0 ;will perform CDS subtraction and extract info from the headers
	sort_bads=1
	initial_badsig=2.
	
do_cuber=0

do_combine=0 ;sorting actually happens here, sorry for poor naming of modules

do_badpix_fix=0
	badsig=2.


do_center=0

do_clean=0
	cross_thresh=0.55

do_inject=1
	;need rplanet in pixels, tplanet in radians, and contrast arrays
	rplanet=[0.3,0.3,0.3]/pixelscale;[0.6]/pixelscale
	tplanet=242.+[90.,0.,180.]+90.
	tplanet=(tplanet+90.)*!DTOR ;convert to N up and radians
	contrast=fltarr(n_elements(rplanet)) & contrast[*]=3.0E-3

	

use_inj=0
start_wave=6 & end_wave=86;82
;start_wave=0 & end_wave=97

comb_type='nwadi'

do_rotate=0
	filter=13.

do_adiklip=0
	do_destripe=1
	default_bin=1 		;auto_bin is turned off below for now, so this default value is used instead
	szz=31.			;half size of processed frame - do not change!!
	k_adiklip=4		;number of KL basis vectors to retain (7)
	wr = 12. 		;Width of annuli in pixels (14) (12 recently)
	nrings = fix(szz/wr) 	;Number of annuli (15)
	n_ang = 1.		;Number of segments (6)
	ANGSEP=0.5	;Exclude frames from ADI KLIP that are within ANGSEP x FWHM of the target frame
	anglemax=360.
	klip_filt=filter ;11.	
	hyp=0
	annmode=1 & if hyp  then annmode=0
	;annmode_inout=[17.-wr+2,17.+wr-2]	;will process only this annulus (in angular segments defined by n_ang)
	annmode_inout=[2,24.]	;will process only this annulus (in angular segments defined by n_ang)

do_sdiklip=1
	k_sdiklip=2.		;number of KL basis vectors to retain (7)
	sdisep=1.5;1.5


suffix=''	

;wavelengths as output by the new cuber
waves=[2.77422155, 2.79857051, 2.82247516, 2.84595897, 2.86904342, 2.89174822,$
 2.91409149, 2.93608999, 2.95775919, 2.97911346, 3.00016614, 3.02092965,$
 3.0414156 , 3.06163483, 3.08159749, 3.10131311, 3.12079064, 3.1400385,$
 3.15906462, 3.17787651, 3.19648124, 3.2148855 , 3.23309564, 3.25111768,$
 3.26895732, 3.28662  ,  3.30411088, 3.32143487, 3.33859667, 3.35560074,$
 3.37245136, 3.38915261, 3.40570839, 3.42212244, 3.43839833, 3.4545395,$
 3.47054923, 3.4864307,  3.50218692, 3.51782083, 3.53333522, 3.54873279,$
 3.56401615, 3.57918779, 3.59425014, 3.60920552, 3.62405617, 3.63880428,$
 3.65345193, 3.66800115, 3.68245389, 3.69681206, 3.71107748, 3.72525192,$
 3.73933712, 3.75333472, 3.76724635, 3.78107358, 3.79481791, 3.80848082,$
 3.82206375, 3.83556807, 3.84899515, 3.86234628, 3.87562275, 3.88882579,$
 3.90195661, 3.91501638, 3.92800624, 3.94092731, 3.95378066, 3.96656734,$
 3.97928839, 3.99194481, 4.00453756, 4.0170676 , 4.02953587, 4.04194326,$
 4.05429065, 4.06657891, 4.07880888, 4.09098139, 4.10309723, 4.11515718,$
 4.12716202, 4.1391125 , 4.15100934, 4.16285326, 4.17464496, 4.18638513,$
 4.19807443, 4.20971353, 4.22130305, 4.23284364, 4.24433591, 4.25578044,$
 4.26717785, 4.2785287 , 4.28983356]


;waves=waves[4:97]


;load pre-determined bad pixel map from Jordan
;badmap=readfits('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/LEECH/median_dark/bpm_and_these_hots.fits')
;badmap[where(badmap eq 0)]=2
;badmap[where(badmap eq 1)]=0
;badmap[where(badmap eq 2)]=1

;dark=readfits('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/LEECH/median_dark/median_cycle_DRK-badfix.fits')
;dark=dark[*,*,1]
;remove_bad_pixels,badmap,'/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/LEECH/median_dark/median_cycle_DRK-badfix.fits'

;sky=readfits('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/LEECH/median_sky/median_cycle_SKY-badfix.fits' )
;remove_bad_pixels,badmap,'/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/LEECH/median_sky/median_cycle_SKY-badfix.fits' 




;------------------------------[ Begin Thought Process  ]---------------------------------


if do_make_skys eq 1 then begin

	files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/raw/','*.fits')

	fpnod=100.

	frames=n_elements(files)

	print, 'Number of frames = ',frames

	print, 'Number of nods = ',frames/fpnod

	
	switchh=0
	z=0
	nn=0
	setflag=1
for ii=0,frames-1 do begin

	print, ii, ' ',files[ii]

	frame=readfits(files[ii],hdr) ;reading in a two-frame file, split into two here
	print, size(frame)
	;if n_elements(size(frame)) gt 5 then begin
	framezero=frame[*,*,0] & frameone=frame[*,*,1]

	flag=strcompress( fxpar(hdr,'FLAG'),/rem)
	print, flag
	
	if flag eq 'NOD_B' then flag='SKY' else flag='PRI'
	;print, 'ii mod fpnod = ', ii mod fpnod 
	;print, 'setflag=',setflag
	;if ii mod fpnod eq 0 and setflag eq 1 then setflag=0 else $
	;if ii mod fpnod eq 0 and setflag eq 0 then setflag=1

	;if setflag then flag='SKY'
	;print, flag

	if flag eq 'SKY' then begin

		if switchh eq 0 then skycubezero=framezero else skycubezero= [ [[skycubezero]], [[framezero]] ]
		if switchh eq 0 then skycubeone=frameone else skycubeone= [ [[skycubeone]], [[frameone]] ]
		switchh = 1
		z=z+1

	endif

	;endif ;frame cube size if
	print, switchh,z,nn
	if switchh eq 1 and z gt 1 and flag ne 'SKY' then begin
		;print, 'Triggered' & hak
		nn=nn+1
		writefits,strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/skys/sky'+string(nn)+'.fits',/rem),[ [[median(skycubezero,dim=3)]],[[median(skycubeone,dim=3)]] ]
		switchh=0
		z=0
	endif

	
endfor






endif ; make sky if


;------------------------------[ Begin Thought Process  ]---------------------------------


if do_sub_sky eq 1 then begin

	files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/raw/','*.fits')

	fpnod=100.

	frames=n_elements(files)

	print, 'Number of frames = ',frames

	print, 'Number of nods = ',frames/fpnod

	
	switchh=0
	z=0
	nn=1
	setflag=1
for ii=0,frames-1 do begin

	print, files[ii]

	frame=readfits(files[ii],hdr)
	framezero=frame[*,*,0] & frameone=frame[*,*,1] ;preserve CDS states
	print, 'Using sky: ',strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/skys/sky'+string(nn)+'.fits',/rem)
	
	
	
	;can improve skysub  with averaging here
	if nn gt 6 then nn=6
	if nn eq 1 or nn ge 6 then sky=readfits(strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/skys/sky'+string(nn)+'.fits',/rem) )
	if nn gt 1 and nn lt 6 then begin
	
		firstsky=readfits(strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/skys/sky'+string(nn-1)+'.fits',/rem) )
	    secondsky=readfits(strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/skys/sky'+string(nn)+'.fits',/rem) )
		fskyzero=firstsky[*,*,0] & fskyone=firstsky[*,*,1]
		sskyzero=secondsky[*,*,0] & sskyone=secondsky[*,*,1]
		
		skyzero=(fskyzero+sskyzero)/2.
		skyone=(fskyone+sskyone)/2.
		

	endif else begin
		skyzero=sky[*,*,0] & skyone=sky[*,*,1]
	endelse
	

	;if nn gt 1 and nn lt 16 then begin ;combine skys when surrounding science frames
	;	sky2=readfits(strcompress('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/skys/sky'+string(nn+1)+'.fits',/rem) )
	;	sky=(sky+sky2)/2.
	;endif

	flag=strcompress( fxpar(hdr,'FLAG'),/rem)


print, flag
	
	if flag eq 'NOD_B' then flag='SKY' else flag='PRI'
;if flag eq 'NOD_B' then flag='SKY'
;	print, 'ii mod fpnod = ', ii mod fpnod 
;	print, 'setflag=',setflag
;	if ii mod fpnod eq 0 and setflag eq 1 then setflag=0 else $
;	if ii mod fpnod eq 0 and setflag eq 0 then setflag=1
;
;	if setflag then flag='SKY' else flag='PRI'
;	print, flag


	if flag eq 'PRI' then begin
		if not no_sky then framezero=framezero-sky
		if not no_sky then frameone=frameone-sky

		;frame=destripe( frame,90.,clip_level=0.0,/nodisp) ;turning off for now, but maybe bring back!

		writefits,strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/skysub/'+obj+'_'+string(ii)+'_'+flag+'.fits',/rem),[ [[framezero]], [[frameone]] ],hdr
		switchh = 1
		;z=z+1
	endif
	print, switchh,z,nn
	if switchh eq 1 and flag ne 'PRI' then begin
		;print, 'Triggered' & hak
		nn=nn+1
		;writefits,strcompress('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/skysub/sky'+string(nn)+'.fits',/rem),median(skycube,dim=3)
		switchh=0
		;z=0
	endif

	
endfor

endif



;------------------------------[ Begin Thought Process  ]---------------------------------


if do_preproc then begin

;first, delete old preproc files
files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/preproc/',count=cnt)
if cnt ge 1 then print, 'Found old preproc files, deleting...'
if cnt ge 1 then file_delete,files,/recursive

;copy sky subtracted data over to new preproc folder
if cnt ge 1 then print, 'Copying sky subtracted data...'
file_copy,'/Volumes/RAID36TB/LMIRCam/TYC5709_2/skysub/','/Volumes/RAID36TB/LMIRCam/TYC5709_2/preproc/',/recursive

;spawn a shell and run the preproc script, which is already pre-set to the correct directory
	 spawn,"python3 '/Volumes/RAID36TB/LMIRCam/TYC5709_2/preproc.py' "


endif ; preproc



if do_sort then begin
	;load primary cubes and combine into single cube

	;cube=fltarr(64,67,9,128) & angles=[]
	files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/preproc/','*.fits')
	angles=[] & flags=[] & dits=[]
	;badmap=readfits('/Volumes/Volumes/RAID36TB/LMIRCam/MWC758_6/bpm_and_these_hots_new3.fits')
	for ii=0,n_elements(files)-1 do begin
		print, 'On file',ii,'/',n_elements(files)-1
		print, files[ii]
		input=readfits(files[ii], hdr)
		

		;neither is needed for preprocessed data (already just takes the final read)
		;input=input[*,*,1] ;- input[*,*,0] ;performs CDS subtraction
		;or...
		;input=input[*,*,1] ;take the final frame from the CDS only 

		;do bad pixel correction here (if uncommented)
		;needs bads = 0 goods = 1
		;input=maskinterp(input,badmap,3,6,"splinterp")


		;end badpix

	
		;writefits,'~/Desktop/test.fits',input
		;if ii eq 0 then cube=input else cube=[ [[cube]], [[input]] ]
		angle=fxpar(hdr,'LBT_PARA')
		dit=fxpar(hdr,'ITIME')
		flag=strcompress(fxpar(hdr,'FLAG'),/rem)
		time=fxpar(hdr,'LBT_UTC')
		print, 'Angle=',angle,' dit= ', dit, 'flag=', flag

		angles=[angles,angle]
		dits=[dits,dit]
		flags=[flags,flag]
		
		
		if sort_bads then begin
		
		if ii eq 0 then begin;only identify bads from first frame to save time
		badmap=input
		print, 'Finding bad pixels'

		badmap[*]=1


		 smallcube = input; make a new bad pixel map for each wavelength
 			cs=2
  			for ix=cs, (size(input))(1)-1-cs do begin
   			 	 for iy=cs, (size(input))(2)-1-cs do begin
   			 	 	print, 'Column ',ix,(size(input))(1)-1-cs
        		   box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs,*]
         		   this =  box[cs,cs,*]  ; this pixel
         		   if mean(this) lt -100 then badmap[ix,iy]=0
         		   box[cs,cs,*] = !values.f_nan
          		   if mean(this) gt (mean(box,/nan)+initial_badsig*stddev(box,/nan)) $
          		   		or mean(this) lt (mean(box,/nan)-initial_badsig*stddev(box,/nan)) $
          		   		then bad=1 else bad=0
          		   if bad then badmap[ix,iy] = 0        
          		   if bad then print, 'Bad pixel found'   
      			  endfor
    	 		  ;print, ' Left Column : ', ix
     		endfor
     		;writefits,'~/Desktop/badmap.fits',badmap
          ;	writefits,'~/Desktop/test.fits',cubefixed[*,*,0,kk]
			print, 'Bad pixels identified = ',float(n_elements(where(badmap eq 0)))
			print, 'Bad pixel fraction = ',float(n_elements(where(badmap eq 0)))/float(n_elements(badmap))

			; BAD PIXEL REMOVAL
  			print
  			;for ii=0, (size(cube))(3)-1 do  begin
				;print, 'Working on frame ', ii, ' / ', (size(cube))(3)-1

		endif ;only find bads on first frame since it takes so long
				writefits,strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/badmap.fits',/rem),badmap,hdr

				print, 'Fixing bad pixels'
  				input = maskinterp(input,badmap,1,3,"splinterp")
  			;endfor
  			;  writefits,'~/Desktop/test-fixed.fits',cubefixed[*,*,0,kk]
	;	hak
		endif
		
		
		
		writefits,strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/cds/'+obj+'_'+string(ii)+'_'+flag+'.fits',/rem),input,hdr
		
		
		;now correct bad pixels
		;remove_bad_pixels,badmap,strcompress('/Volumes/RAID36TB/LMIRCam/TYC5709_2/cds/'+obj+'_'+string(ii)+'_'+flag+'.fits',/rem)
	
		;now subtract the median dark or sky (or none)
		
		;cds=readfits(strcompress('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/cds/'+obj+'_'+string(time)+'_'+flag+'.fits',/rem),hdr)
		
		;cds=cds-dark
		;sky=sky-dark
		;cds=cds-sky
		;writefits,strcompress('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/cds/'+obj+'_'+string(time)+'_'+flag+'.fits',/rem),cds,hdr

		;now do CDS
		;cds=readfits(strcompress('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/cds/'+obj+'_'+string(ii)+'.fits',/rem),hdr)
		;cds=cds[*,*,1]-cds[*,*,0]
		;writefits,strcompress('/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/cds/'+obj+'_'+string(ii)+'.fits',/rem),cds,hdr
	endfor

	;cube=cube[0:63,0:63,*,26:125]

	;angles=[-36.7287, -31.762, -26.09, -19.74, -11.98, -4.52, 3.1, 10.62, 28.817]

	;writefits,strcompress(folder+obj+'_cube.fits',/rem),cube
	;save, filename=strcompress(folder+obj+'_angles.sav',/rem),angles,dits,flags

	;print, angles

endif ;do_sort if



;------------------------------[ Begin Thought Process  ]---------------------------------


if do_cuber then begin

;files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/incubes/',count=cnt)
;if cnt ge 1 then print, 'Found old incubes files, deleting...'
;if cnt ge 1 then file_delete,files,/recursive

;copy sky subtracted data over to new preproc folder
;if cnt ge 1 then print, 'Copying sorted data...'
;file_copy,'/Volumes/RAID36TB/LMIRCam/TYC5709_2/cds/','/Volumes/RAID36TB/LMIRCam/TYC5709_2/incubes/',/recursive


		files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/incubes/','*.fits',count=cnt)
		if cnt gt 1 then file_delete,files

		spawn,"python3 '/Users/kevinwagner/idl/cuber_wrapper_wispit2.py' "
endif

;------------------------------[ Begin Thought Process  ]---------------------------------

if do_combine then begin

	files=file_search('/Volumes/RAID36TB/LMIRCam/TYC5709_2/incubes/','*.fits')

	cube=fltarr(67,63,n_elements(files),99)

	starti=0

	for ii=starti, n_elements(files)-1 do begin
		print, ii,'/',n_elements(files)-1
		c=readfits(files[ii],hdr)
		cube[*,*,ii,*]=c
		
		angle=float(fxpar(hdr,'LBT_PARA'))
		print, angle
		;hak
		exptime=fxpar(hdr,'EXPTIME')
		flag=strcompress(fxpar(hdr,'FLAG'),/rem)
		RA=fxpar(hdr,'LBT_RA')
		DEC=fxpar(hdr,'LBT_DEC')


		if ii eq starti then angles=angle else angles=[angles,angle]
		if ii eq starti then exptimes=exptime else exptimes=[exptimes,exptime]
		if ii eq starti then flags=flag else flags=[flags,flag]
		if ii eq starti then RAs=RA else RAs=[RAs,RA]
		if ii eq starti then DECs=DEC else DECs=[DECs,DEC]
		;print, angles, exptimes, flags, RAs, DECs
	endfor

	;cube=cube[*,3:66,*,42:116]
	;cube=cube[*,*,*,4:97]

	pricube=cube
	;pricube=cube[*,*,where(flags eq 'PRI'),*]
	;skycube=cube[*,*,where(flags eq 'SKY'),*]
	;darkcube=cube[*,*,where(flags eq 'DRK'),*]

	;darkmed=median(darkcube,dim=3)
	;skymed=median(skycube,dim=3)
	;writefits,strcompress(folder+obj+'_dark.fits',/rem),darkmed
	;writefits,strcompress(folder+obj+'_sky.fits',/rem),skymed

	;sorting by parallactic angle
	
	
	if 0 then begin
	order=sort(angles)
	angles=angles[order]
	exptimes=exptimes[order]
	flags=flags[order]
	RAs=RAs[order]
	DECs=DECs[order]
	pricube=pricube[*,*,order,*]
	endif

	print,angles
	hak

	writefits,strcompress(folder+obj+'_pricube.fits',/rem),pricube
	;writefits,strcompress(folder+obj+'_skycube.fits',/rem),skycube
	;writefits,strcompress(folder+obj+'_darkcube.fits',/rem),darkcube
	save,filename=strcompress(folder+obj+'_info.sav',/rem),angles,exptimes,flags,RAs,DECs

	;angles=angles[where(flags eq 'PRI')]

	save,filename=strcompress(folder+obj+'_primary_angles.sav',/rem),angles

endif

if do_badpix_fix then begin

	cube=readfits(strcompress(folder+obj+'_pricube.fits',/rem))
	cubefixed=cube
	for kk=0,(size(cube))(4)-1 do begin ;loop through waves
		
		print, 'Finding bad pixels for wave',kk
		badmap=cube[*,*,0,kk]

		badmap[*]=1


		 smallcube = cube[*,*,0:50,kk]; make a new bad pixel map for each wavelength
 			cs=2
  			for ix=cs, (size(cube))(1)-1-cs do begin
   			 	 for iy=cs, (size(cube))(2)-1-cs do begin
        		   box = smallcube[ix-cs:ix+cs,iy-cs:iy+cs,*]
         		   this =  box[cs,cs,*]  ; this pixel
         		   if mean(this) lt -100 then badmap[ix,iy]=0
         		   box[cs,cs,*] = !values.f_nan
          		   if mean(this) gt (mean(box,/nan)+badsig*stddev(box,/nan)) $
          		   		or mean(this) lt (mean(box,/nan)-badsig*stddev(box,/nan)) $
          		   		then bad=1 else bad=0
          		   if bad then badmap[ix,iy] = 0           
      			  endfor
    	 		  ;print, ' Left Column : ', ix
     		endfor
     		;writefits,'~/Desktop/badmap.fits',badmap
          ;	writefits,'~/Desktop/test.fits',cubefixed[*,*,0,kk]


			; BAD PIXEL REMOVAL
  			print
  			for ii=0, (size(cube))(3)-1 do  begin
				;print, 'Working on frame ', ii, ' / ', (size(cube))(3)-1
  				cubefixed[*,*,ii,kk] = maskinterp(cube[*,*,ii,kk],badmap,1,3,"splinterp")
  			endfor
  			;  writefits,'~/Desktop/test-fixed.fits',cubefixed[*,*,0,kk]
	;	hak
		med=median(cubefixed[*,*,*,kk],dimension=3)
		if kk eq 0 then medcube=med else medcube = [ [[medcube]], [[med]] ]
	endfor

	writefits,strcompress(folder+obj+'_pupil_badfix_cube.fits',/rem),medcube
  		
  	writefits,strcompress(folder+obj+'_pricube_badfix.fits',/rem),cubefixed


endif

if do_center then begin

	cube=readfits(strcompress(folder+obj+'_pricube_badfix.fits',/rem))
	cube[where(finite(cube) eq 0)]=0.
	half_sizex=31. & half_sizey=31.
	cube=cube[0:61,0:61,*,*]
	cencube=cube
	
	for kk=0,(size(cube))(4)-1 do begin ;loop through waves
		for ii=0,(size(cube))(3)-1 do begin

			cencube[*,*,ii,kk]=cencube[*,*,ii,kk]-smooth(cencube[*,*,ii,kk],9.)
			cencube[*,*,ii,kk]=smooth(cencube[*,*,ii,kk],5.)
			
		endfor
	endfor
	
	
	for kk=0,(size(cube))(4)-1 do begin ;loop through waves

		first=cencube[*,*,0,kk] ;pick the first frame
		;peak=mpfit2dpeak(first,A,/tilt) ;fit the peak
		;xx=A[4] & yy=A[5]	;extract coords
		cntrd,first,33,33,xx,yy,2
		cencube[*,*,0,kk]=fshift(first,30.5-xx,30.5-yy) 
		cube[*,*,0,kk]=fshift(cube[*,*,0,kk],30.5-xx,30.5-yy) 

		print, xx,yy

		for ii=0,(size(cube))(3)-1 do begin
			
			corr=crosscorr(cencube[*,*,0,kk],cencube[*,*,ii,kk],pmax, /cut)
			if pmax(0) ge 0. then dx = half_sizex-pmax(0) else dx = -half_sizex+abs(pmax(0))
			if pmax(1) ge 0. then dy = half_sizey-pmax(1) else dy = -half_sizey+abs(pmax(1))

			print, 'Shifting frame', ii, ' / ',(size(cube))(3)-1,' of kk=',kk,' by: ',String(-dx),String(-dy)
			cube[*,*,ii,kk]=fshift(cube[*,*,ii,kk],-dx,-dy)
		endfor

		med=median(cube[*,*,*,kk],dimension=3)
		;writefits,'~/Desktop/test.fits',med
		if kk eq 0 then medcube=med else medcube = [ [[medcube]], [[med]] ]
	endfor

	writefits,strcompress(folder+obj+'_cube_cen.fits',/rem),cube
	writefits,strcompress(folder+obj+'_pupil_cube.fits',/rem),medcube

endif ;center if

if do_clean then begin

	cube=readfits(strcompress(folder+obj+'_cube_cen.fits',/rem))
	restore,filename=strcompress(folder+obj+'_primary_angles.sav',/rem);,angles


;performing auto bad frame recognition here
;first calculate median of science cube in wavelength
for ii=0,n_elements(angles)-1 do begin
	medframe=median(reform(cube[*,*,ii,*]),dim=3)
	;if pup_stab_clean then medframe=rot(medframe,-pa[ii]-pupoff[ii],/interp)
	if ii eq 0 then medframes=medframe else medframes=[ [[medframes]],[[medframe]] ]

endfor


;then compute median of those frames

medarr, medframes, pup_median

;create gaussian instead
pup_median=medframes[*,*,200]
pup_median=gauss2dfit(pup_median)


sz=(size(pup_median))(1) 
width=4.
pup_median = psf_Gaussian(NPIX=sz, FWHM=[width,width])

writefits,folder+obj+'_cenref.fits',pup_median

;block center

for xx=0,61 do for yy=0,61 do if sqrt( (xx-31.)^2. + ( (yy-31.)^2.)) lt 5 then pup_median[xx,yy]=0.

;now check frames against the median
for ii=0,n_elements(angles)-1 do begin
	medframe=median(reform(cube[*,*,ii,*]),dim=3)
	for xx=0,61 do for yy=0,61 do if sqrt( (xx-31.)^2. + ( (yy-31.)^2.)) lt 5 then medframe[xx,yy]=0.
	;if pup_stab_clean then medframe=rot(medframe,-pa[ii]-pupoff[ii],/interp)
	maxcor=max(crosscorr(medframe,pup_median))
			print, 'Frame ',ii,' has max correlation value of ', maxcor, ' with pupil median.'
		if ii eq 0 then corrs=maxcor else corrs=[corrs,maxcor]
endfor

;keep those above cross_thresh percentile (defined out of unity), e.g. 0.8=80%

	;lcube=lcube[*,*,where(corrs ge cross_thresh)]
	;rcube=rcube[*,*,where(corrs ge cross_thresh)]
;cube=cube[*,*,*,where(corrs ge cross_thresh)]
	

	print, 'Cleaned ', n_elements(where(corrs lt cross_thresh)) ,' frames.'
	;print, 'Bad frames:', where(corrs lt cross_thresh)

bads=where(corrs lt cross_thresh)

if  n_elements(where(corrs lt cross_thresh)) ge 0.5*n_elements(corrs) then begin
		print, '!!! Warning! Autoclean has rejected over 50% of frames !!! 
		print, 'You will porbably want to go back to the reduce_ifs.pro routine and change the cross_thresh parameter to a lower value (0 will accept all frames).'
		print, 'If this is a binary star, this may have happened because of comparing single frames (with a bright star) to a median-combined pupil containing a smeared image of the star. Try enabling pup_stab_clean in this case.'
		print, 'It is recommended to apply the /skip keyword next run to reduce tiem spet on initial steps.'
		print, 'Hit any key to continue anyways, or ESC to stop the reduction.'
		hak
	endif

goods=fltarr((size(angles))(1))
goods[*]=1
goods[bads]=0

goods[where(corrs lt cross_thresh)]=0

goods=where(goods gt 0)

cube=cube[*,*,goods,*]

angles=angles[goods]

	writefits,strcompress(folder+obj+'_cube_cen_clean.fits',/rem),cube

	save,filename=strcompress(folder+obj+'_primary_angles_clean.sav',/rem),angles


endif ;clean if



if do_inject  then begin


;need rplanet in pixels, tplanet in radians, and contrast arrays

scicube=readfits(strcompress(folder+obj+'_cube_cen_clean.fits',/rem),hdr)
restore,strcompress(folder+obj+'_primary_angles_clean.sav',/rem)

szz=(size(scicube))(1)/2. & cnt=(size(scicube))(3) & nwaves=n_elements(waves) & nplanets=n_elements(contrast)

bigarr=fltarr(4.*szz,4.*szz,cnt,nwaves)
bigpsfarr=fltarr(4.*szz,4.*szz,nwaves)
psf=readfits('/Volumes/RAID36TB/LMIRCam/TYC5709_2/reduced/PDS201_xyl_derot-nf.fits')
print, size(bigpsfarr[szz:3.*szz-1,szz:3.*szz-1,*])
print, size(psf)
print, size(bigpsfarr)
bigpsfarr[szz:3.*szz-1,szz:3.*szz-1,*]=reform(psf)
print, szz,3.*szz-1,szz,3.*szz-1
for k=0,cnt-1 do begin


	for ii=0, nwaves-1 do begin	
		
			if k eq 0 then gpsf=gauss2dfit(bigpsfarr[*,*,ii],/tilt)
			if k eq 0 then opsf=bigpsfarr[*,*,ii]
			;writefits,'~/Desktop/gpsf.fits',[[[opsf]],[[gpsf]] ]


			bigarr[szz:3.*szz-1,szz:3.*szz-1,k,ii]=scicube[*,*,k,ii]

			bigarr[*,*,k,ii]=rot(bigarr[*,*,k,ii],0.,2.0,/interp)  ;blow up
 			bigarr[*,*,k,ii]=rot(bigarr[*,*,k,ii],-angles[k],cubic=-0.5);/interp),/INTERP)
		 				
			spsf=reform(bigarr[*,*,k,ii]) ;use same frame
 				
 					xcc=31.
 					ycc=31.
 					if k eq 0 and ii eq 0 then print, $
 						'-------- Adding synthetic planets -------------'
 			
 				for jj=0,nplanets-1 do begin 

				print, 'Injecting planet # ', jj,ii,k, ' out of ', nplanets-1, nwaves-1, cnt-1				
 			 				
 				bigarr[*,*,k,ii]= bigarr[*,*,k,ii] + $
 					fshift(contrast[jj]*spsf, -rplanet[jj]*2.0*cos(tplanet[jj]),-rplanet[jj]*2.0*sin(tplanet[jj]))

				 print,'Injection coordinates: ',31.-rplanet[jj]*cos(tplanet[jj]),31.-rplanet[jj]*sin(tplanet[jj])
 				
 				endfor
				if ii eq 0 then print, 'Please wait...'
 			 	bigarr[*,*,k,ii]= rot(bigarr[*,*,k,ii],angles[k],cubic=-0.5);/interp),/INTERP)
				bigarr[szz:3.*szz-1,szz:3.*szz-1,k,ii]=congrid(bigarr[*,*,k,ii],2.*szz,2.*szz) ;shrink down
 			 	scicube[*,*,k,ii]=bigarr[szz:3.*szz-1,szz:3.*szz-1,k,ii]
				
		if k eq cnt-1 then begin
			med=median(scicube[*,*,*,ii],dim=3)
			if ii eq 0 then meds=med else meds=[[[meds]],[[med]]]
		endif
	endfor
endfor 	

writefits,strcompress(folder+obj+'_cube_cen_clean_inj.fits',/rem),scicube,hdr
writefits,strcompress(folder+obj+'_pupil_cen_clean_inj.fits',/rem),meds,hdr
endif ;addplanets if



if do_rotate then begin

	offset=180.

	cube=readfits(strcompress(folder+obj+'_cube_cen_clean.fits',/rem))

	if use_inj then cube=readfits(strcompress(folder+obj+'_cube_cen_clean_inj.fits',/rem))

	restore,strcompress(folder+obj+'_primary_angles_clean.sav',/rem)
	


nwaves=n_elements(waves)
cube=cube[*,*,*,start_wave:end_wave]
waves=waves[start_wave:end_wave]
nwaves=n_elements(waves)


	if do_destripe then begin

		for i=0,(size(cube))(3)-1 do for j=0,nwaves-1 do cube[*,*,i,j]=destripe(cube[*,*,i,j],90.,clip_level=0.0,/nodisp)
		for i=0,(size(cube))(3)-1 do for j=0,nwaves-1 do cube[*,*,i,j]=destripe(cube[*,*,i,j],0.,clip_level=0.0,/nodisp)

	endif
	if filter gt 0 then for kk=0,(size(cube))(4)-1 do for ii=0,(size(cube))(3)-1 do cube[*,*,ii,kk]=cube[*,*,ii,kk]-smooth(cube[*,*,ii,kk],filter)


	 subcube=cube & rotcube=cube






	for kk=0,(size(cube))(4)-1 do begin ;loop through waves
		;compute a median
		;med=median(cube[*,*,*,kk],dimension=3)
		medarr,cube[*,*,*,kk],med
		;writefits,'~/Desktop/test.fits',med
	for ii=0,(size(cube))(3)-1 do begin	;loop through frames
		;subract it
		subcube[*,*,ii,kk]=reform(cube[*,*,ii,kk])-reform(med)
		;derotate
		rotcube[*,*,ii,kk]=rot(cube[*,*,ii,kk],-angles[ii]-offset,/interp)
		subcube[*,*,ii,kk]=rot(subcube[*,*,ii,kk],-angles[ii]-offset,/interp)
	endfor

		;make a median
		if comb_type eq 'mean' then adimed=mean(subcube[*,*,*,kk],dimension=3)
		if comb_type eq 'median' then adimed=median(subcube[*,*,*,kk],dimension=3)
		if comb_type eq 'nwadi' then adimed=nw_ang_comb(subcube[*,*,*,kk],angles)
		rotmed=median(rotcube[*,*,*,kk],dimension=3)

		;store in cubes
		if kk eq 0 then adimedcube=adimed else adimedcube = [ [[adimedcube]], [[adimed]] ]
		if kk eq 0 then rotmedcube=rotmed else rotmedcube = [ [[rotmedcube]], [[rotmed]] ]

	endfor


	writefits,strcompress(folder+obj+'_xyln_adi.fits',/rem),subcube
	writefits,strcompress(folder+obj+'_xyln_derot.fits',/rem),rotcube


	writefits,strcompress(folder+obj+'_xyl_adi.fits',/rem),adimedcube
	writefits,strcompress(folder+obj+'_xyl_derot.fits',/rem),rotmedcube

	adimed=mean(adimedcube,dimension=3)
	rotmed=mean(rotmedcube,dimension=3)

	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) lt annmode_inout[0] then adimed[xx,yy]=0.
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) gt annmode_inout[1] then adimed[xx,yy]=0.
	writefits,strcompress(folder+obj+'_adi.fits',/rem),adimed
find_sources,strcompress(folder+obj+'_adi'+suffix+'.fits',/rem),fwhm=(3.5*1E-6) / (8.4) * 206265. / 0.0345
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) lt annmode_inout[0] then rotmed[xx,yy]=0.
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) gt annmode_inout[1] then rotmed[xx,yy]=0.
	writefits,strcompress(folder+obj+'_derot.fits',/rem),rotmed
find_sources,strcompress(folder+obj+'_derot'+suffix+'.fits',/rem),fwhm=(3.5*1E-6) / (8.4) * 206265. / 0.0345

endif ;rotate if



;---------------[ Beginning Reduction Script ]-------------------------------------

if do_adiklip then begin


cube=readfits(strcompress(folder+obj+'_cube_cen_clean.fits',/rem))
	if use_inj then cube=readfits(strcompress(folder+obj+'_cube_cen_clean_inj.fits',/rem))

	restore,strcompress(folder+obj+'_primary_angles_clean.sav',/rem)
;print, angles
;hak

nwaves=n_elements(waves)
cube=cube[*,*,*,start_wave:end_wave]
if nwaves gt end_wave-start_wave+1 then waves=waves[start_wave:end_wave]
nwaves=n_elements(waves)
cnt=n_elements(angles)
if do_destripe then begin

for i=0,cnt-1 do for j=0,nwaves-1 do cube[*,*,i,j]=destripe(cube[*,*,i,j],90.,clip_level=0.0,/nodisp)
for i=0,cnt-1 do for j=0,nwaves-1 do cube[*,*,i,j]=destripe(cube[*,*,i,j],0.,clip_level=0.0,/nodisp)

endif



;bins to about a size of 20

auto_bin=0 & nwaves=n_elements(waves)

cube=cube[*,*,sort(angles),*]
angles=angles[sort(angles)]

if auto_bin then begin
	bin=1
	if n_elements(angles) gt 40 then bin=2.
	if n_elements(angles) gt 80 then bin=4.
	if n_elements(angles) gt 120 then bin=6.
	if n_elements(angles) gt 160 then bin=8.
	if n_elements(angles) gt 200 then bin=10.
endif else bin = default_bin

ncubesize=fix(float(n_elements(angles))/float(bin))
ncube=fltarr(62,62,ncubesize,nwaves)
print, 'New cube size will be: ',size(ncube)


precnt=n_elements(angles)
if bin gt 1 then begin
	for kk=0,nwaves-1 do begin
	binnum=0
	;st=1
	for ii=0.,(size(cube))(3)-1. do begin
		;print, fix(ii+1.) mod fix(bin)
		if fix(ii+1.) mod fix(bin) eq 1 then binned=cube[*,*,ii,kk] else binned=[ [[binned]], [[ reform(cube[*,*,ii,kk]) ]] ]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
		

			 binned=mean(binned,dim=3)
			

			;print, kk, binnum
			ncube[*,*,binnum,kk]=binned
			binnum=binnum+1
		endif
	endfor	
	;print, size(binned_cube)
	;kcube=binned_cube
	endfor

		bin=float(bin)


	;st=1
	binned_angles=[]
	for ii=0,precnt-1. do begin
		print, angles[ii]
		if fix(ii+1.) mod fix(bin) eq 1 then binned_angle=angles[ii] else binned_angle=[binned_angle,float(angles[ii])]
		if fix(ii+1.) mod fix(bin) eq 0 then begin
			;if st eq 1 then begin 
					binned_angles = [binned_angles,mean(binned_angle)]; binned_angle/bin & st=0 & 
			;endif else binned_angles=[binned_angles,binned_angle/bin]			
		endif
	endfor	
	angles=binned_angles
cube=ncube
endif

;print, angles
;hak

;end binning code


offset=180.
;angles=angles+offset


final_cube=cube
	save,filename=strcompress(folder+obj+'_primary_angles_clean_binned.sav',/rem),angles
sangles=angles



;waves=fltarr(nwaves) & waves[*]=4;E-6

cnt=(size(cube))(3)

if klip_filt gt 1 then for ii=0,cnt-1 do for jj=0,nwaves-1 do cube[*,*,ii,jj]=cube[*,*,ii,jj]-smooth(cube[*,*,ii,jj],klip_filt)


for kk=0, nwaves-1 do begin

kcube=reform(cube[*,*,*,kk])
;print,n_elements(where(finite(kcube) lt 1))

kcube[where(finite(kcube) lt 1)]=0. ; replace NaNs
writefits,strcompress(folder+obj+'_klip_input_cube.fits'),kcube,hdr
klip=kcube
;this loops through the frames in each particular band and performs the KLIP routine
	for ii=0, cnt-1 do begin
		angles=sangles
		print, 'input size:', size(kcube)
		print, '--[ KLIPing Image ', ii+1, ' / ', cnt,']--'
		print, 'On wavelength:', kk+1, ' / ', nwaves, '     =     ',waves[kk],' µm '
			
  		if hyp  then  klip = adiklip(kcube, k_adiklip, target=ii, $
  		 anglemask=anglemask, distmask=distmask, posang=-angles,  $
 		 	wl=waves[kk],diam=8.4,pixelscale=0.035,angsep=angsep,anglemax=anglemax,obj=obj,nrings=nrings,wr=wr,$
  		  n_ang =n_ang,spot_radius=spot_radius,rho=rho,phi=phi,/hyper) 
  		 
  		if annmode  then klip = adiklip(kcube, k_adiklip, target=ii, anglemask=anglemask, distmask=distmask,$
  		  posang=-angles,  wl=waves[kk],diam=8.4,pixelscale=0.035,angsep=angsep,anglemax=anglemax,$
  		  obj=obj,nrings=nrings,wr =wr, n_ang =n_ang,annmode_inout=annmode_inout)
  		  
  		  if annmode eq 0 and hyp eq 0 then klip = adiklip(kcube, k_adiklip, target=ii, anglemask=anglemask, 	$
  		  distmask=distmask, posang=-angles,  wl=waves[kk],diam=8.4,pixelscale=0.035,angsep=angsep,anglemax=anglemax,$
  		  obj=obj,nrings=nrings,wr =wr, n_ang =n_ang)
  		 	print, 'output size:', size(klip)

		if ii eq 0 then klipcube=klip else klipcube=[[[klipcube]],[[klip]]]
	endfor
angles=sangles
writefits,strcompress(folder+obj+'_klip_kcube.fits'),klipcube,hdr
	
	final_cube[*,*,*,kk]=klipcube
	
	;save preliminary output
	
		for ii=0, cnt-1 do begin
 			klipcube[*,*,ii]=rot(klipcube[*,*,ii],-angles[ii]-offset,/interp)	
		endfor
		
		
		if comb_type eq 'median' then medarr, klipcube, kframemed
		if comb_type eq 'nwadi' then kframemed=nw_ang_comb(klipcube,angles)
		if comb_type eq 'mean' then kframemed=mean(klipcube,dim=3)

		
		
		writefits,strcompress(folder+obj+'_'+string(kk)+'_prelim.fits',/rem),kframemed,hdr

		if kk eq 0 then medcube=kframemed else medcube=[ [[medcube]], [[kframemed]] ]

endfor

	writefits,strcompress(folder+obj+'_adiklip_xyln.fits'),final_cube,hdr

	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) lt annmode_inout[0] then medcube[xx,yy,*]=!values.f_nan
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) gt annmode_inout[1] then medcube[xx,yy,*]=!values.f_nan
	writefits,strcompress(folder+obj+'_adiklip_xyl.fits'),medcube,hdr
	;medarr,medcube,medframe
	
	
	
	adiklip=mean(medcube,dimension=3)
	adiklipl=mean(medcube[*,*,where(waves gt 3.6 and waves lt 4.0)],dimension=3)
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) lt annmode_inout[0] then adiklip[xx,yy]=!values.f_nan
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) gt annmode_inout[1] then adiklip[xx,yy]=!values.f_nan

	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) lt annmode_inout[0] then adiklipl[xx,yy]=!values.f_nan
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) gt annmode_inout[1] then adiklipl[xx,yy]=!values.f_nan
		writefits,strcompress(folder+obj+'_adiklip.fits'),adiklip,hdr
		writefits,strcompress(folder+obj+'_adiklip_L.fits'),adiklipl,hdr
		
		find_sources,strcompress(folder+obj+'_adiklip'+suffix+'.fits',/rem),fwhm=(3.5*1E-6) / (8.4) * 206265. / 0.0345
		;make a convolved cube
conv_xylcube=medcube
for ccc=0,nwaves-1 do begin
	
	sz=62.
	width=(waves[ccc]*1E-6) / (8.4) * 206265. / 0.035
	print, 'PSF Width: ',width
	PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	PSFN = PSF/MAX(PSF)
	conv_xylcube[where(finite(conv_xylcube) ne 1)]=0.
	conv_xylcube[*,*,ccc] = convolve(conv_xylcube[*,*,ccc], PSFN)

endfor

	writefits,strcompress(folder+obj+'_adiklip_xyl_conv.fits'), conv_xylcube, hdr
	
	;medarr, conv_xylcube, convmed
	
			adiklip=mean(conv_xylcube,dimension=3)

	
		writefits,strcompress(folder+obj+'_adiklip_conv.fits'), adiklip, hdr


oldfiles=FILE_SEARCH(folder,'*prelim.fits',COUNT=cc)
FILE_DELETE, oldfiles




endif ;adiklip


;---------------------------------------------------------------

if do_sdiklip  then begin

debug=0
fresh_sdi=0
comb_type='nw-mean'
suffix=''

if fresh_sdi  then adicube=scicube else $
	adicube=readfits(folder+obj+'_adiklip_xyln.fits')
	
	restore,strcompress(folder+obj+'_primary_angles_clean_binned.sav',/rem)


;adicube=adicube[*,*,start_wave:end_wave,*]
nwaves=n_elements(waves)
if nwaves gt end_wave-start_wave+1 then waves=waves[start_wave:end_wave]


	nwaves=n_elements(waves)
	offset=180.
	cnt=n_elements(angles)


if fresh_sdi  then begin
for k=0,nwaves do begin
		for ii=0, cnt-1 do begin
 		adicube[*,*,ii,k]=rot(adicube[*,*,ii,k],-angles[ii]-offset,/INTERP)
		endfor
endfor
endif


finalcube=adicube ;to be used later
for i=0,cnt-1 do begin

	tempcube=reform(adicube[*,*,i,*]) ;this is a x,y,l cube from just a single IFS frame
	tempcubesc=tempcube	;this will be the scaled cube
	
	print, size(tempcubesc)
	scaling=waves
	
	if debug  then writefits, '~/Desktop/sdi_input_prescaling.fits',tempcube, head
	
	for k=0,nwaves-1 do begin ;now for each lambda slice we scale the cube then perform SDI
		for kk=0,nwaves-1 do begin ;scaling the cube
			tempcubesc[*,*,kk]=rot(tempcube[*,*,kk],0., (scaling[k]/scaling[kk]), /INTERP)
		endfor
		
		tempcubesc[where(finite(tempcubesc) eq 0)]=0.
		starttime=double(systime(/JULIAN))

		if debug  then writefits, '~/Desktop/sdi_input_debug.fits',tempcubesc, head
	
		print, '--[ SDI-KLIPing Image ', i, ' out of ',cnt-1,' ]--'
		print, 'Object: ',obj
		print, 'On wavelength:', waves[k],'/',waves[nwaves-1]
		;now use SDI klip
		if annmode  then tempcubesck = sdiklip(tempcubesc, k_sdiklip, target=k, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=waves[k], diam=8.4, pixelscale=0.0345, angsep=sdisep, obj=obj, nrings=nrings,wr =wr, n_ang =n_ang,annmode_inout=annmode_inout)
		if hyp  then tempcubesck = sdiklip(tempcubesc, k_sdiklip, target=k, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=waves[k], diam=8.4, pixelscale=0.0345, angsep=sdisep, obj=obj, nrings=nrings,wr =wr, n_ang =n_ang,spot_radius=spot_radius,rho=rho,phi=phi, /hyper)
		if annmode eq 0 and hyp eq 0 then tempcubesck = sdiklip(tempcubesc, k_sdiklip, target=k, anglemask=anglemask, distmask=distmask,posang=-angles, scaling=scaling,  wl=waves[k], diam=8.4, pixelscale=0.0345, angsep=sdisep, obj=obj, nrings=nrings,wr =wr, n_ang =n_ang)
		print, '------------------------[ Percent complete:', ((float(i)/float(cnt-1.))+float(k+1.)/(nwaves*float(cnt)))*100.
		endtime=double(systime(/JULIAN))
		print, '------------------------[ Time remaining (min):', (1.-((float(i)/float(cnt-1.))+float(k+1.)/(nwaves*float(cnt))))*(nwaves*float(cnt))*(endtime-starttime) *86400./60.
		
				if debug  then	writefits, '~/Desktop/sdi_output_debug.fits',tempcubesck, head

		finalcube[*,*,i,k]=tempcubesck
		finalcubederot=finalcube
		if k eq nwaves-1 then begin
			;medarr, reform(finalcubederot[*,*,i,*]), framei
			framei=mean(reform(finalcubederot[*,*,i,*]),dim=3)
			framei=rot(framei,-angles[i]-offset,/interp)
			writefits, strcompress(folder+''+obj+'_sdiklip_'+String(i)+'_prelim.fits',/rem),framei, head
			if i eq 0 then sdis=framei else sids=[[[sdis]],[[framei]]]
		endif
	endfor

	

endfor

	finalcubederot=finalcube
	

		;derotates the frames and median combines for each band
	for ii=0,cnt-1 do for k=0,nwaves-1 do finalcubederot[*,*,ii,k]=rot(finalcubederot[*,*,ii,k],-angles[ii]-offset,/interp)
	for k=0,nwaves-1 do begin
		if comb_type eq 'median' then medarr, reform(finalcubederot[*,*,*,k]), kframemed
		if comb_type eq 'mean' then kframemed=mean(reform(finalcubederot[*,*,*,k]),dim=3)
		if comb_type eq 'nw-mean' then kframemed=nw_ang_comb(reform(finalcubederot[*,*,*,k]),angles+offset)
		if k eq 0 then xylcube=kframemed else xylcube=[ [[xylcube]] , [[kframemed]] ]
		if k eq 0 then sumimg=kframemed else sumimg=sumimg+kframemed
	endfor
	
	simimg=sumimg/nwaves
	
	;medarr, xylcube, medimg
	medimg=mean(xylcube,dim=3)
	imgcube=[]
for kz=0, (size(finalcube))(4)-1 do begin &	kzcube=reform(finalcube[*,*,*,kz]) & medarr,kzcube,kzmed & imgcube=[[[imgcube]],[[kzmed]]] & endfor	

;final high pass filter
;xylcube[where(finite(xylcube) eq 0)]=0.
;for k=0,nwaves-1 do xylcube[*,*,k]=xylcube[*,*,k]-smooth(xylcube[*,*,k],filter)

;make a convolved cube
	conv_xylcube=xylcube
imgsdi=mean(xylcube,dimension=3)
writefits, folder+obj+'_sdiklip'+suffix+'.fits', imgsdi, head

new=imgsdi

xylcubel=xylcube[*,*,where(waves gt 3.6 and waves lt 4.)]
imgsdil=mean(xylcubel,dimension=3)
writefits, folder+obj+'_sdiklip_L'+suffix+'.fits', imgsdil, head



;xylcubei=xylcube[*,*,where(waves gt 2.9 and waves lt 3.0)]
;imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_2.9-3.0'+suffix+'.fits', imgsdii, head


;xylcubei=xylcube[*,*,where(waves gt 3.0 and waves lt 3.1)]
;imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.0-3.1'+suffix+'.fits', imgsdii, head




;xylcubei=xylcube[*,*,where(waves gt 3.1 and waves lt 3.2)]
;imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.1-3.2'+suffix+'.fits', imgsdii, head
;bcube=imgsdii
;bcube=[[[bcube]],[[imgsdii]]]

xylcubei=xylcube[*,*,where(waves gt 3.2 and waves lt 3.3)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.2-3.3'+suffix+'.fits', imgsdii, head
bcube=imgsdii
;bcube=[[[bcube]],[[imgsdii]]]

xylcubei=xylcube[*,*,where(waves gt 3.3 and waves lt 3.4)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.3-3.4'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 3.4 and waves lt 3.5)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.4-3.5'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 3.5 and waves lt 3.6)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.5-3.6'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 3.6 and waves lt 3.7)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.6-3.7'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 3.7 and waves lt 3.8)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.7-3.8'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 3.8 and waves lt 3.9)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.8-3.9'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 3.9 and waves lt 4.0)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_3.9-4.0'+suffix+'.fits', imgsdii, head

bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 4.0 and waves lt 4.1)]
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_4.0-4.1'+suffix+'.fits', imgsdii, head
bcube=[[[bcube]],[[imgsdii]]]
xylcubei=xylcube[*,*,where(waves gt 4.1 and waves lt 4.2)]
print,waves
imgsdii=mean(xylcubei,dimension=3)
;writefits, folder+obj+'_sdiklip_4.1-4.2'+suffix+'.fits', imgsdii, head
bcube=[[[bcube]],[[imgsdii]]]
writefits, folder+obj+'_sdiklip_0.1um_cube'+suffix+'.fits', bcube, head
find_sources,strcompress(folder+obj+'_sdiklip'+suffix+'.fits',/rem),fwhm=(3.5*1E-6) / (8.4) * 206265. / 0.0345




for ccc=0,nwaves-1 do begin
	
	sz=62.-1.
	width=(waves[ccc]*1E-6) / (8.4) * 206265. / 0.0345
	print, 'PSF Width: ',width
	PSF = psf_Gaussian(NPIX=sz, FWHM=[width,width])
	PSFN = PSF/MAX(PSF)
	conv_xylcube[where(finite(conv_xylcube) ne 1)]=0.
	conv_xylcube[*,*,ccc] = convolve(conv_xylcube[*,*,ccc], PSFN)

endfor
	
writefits, folder+obj+'_sdiklip_xyln.fits', finalcube, head
writefits, folder+obj+'_sdiklip_xyl'+suffix+'.fits', xylcube, head

imgsdi=mean(conv_xylcube,dimension=3)

for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) lt annmode_inout[0] then imgsdi[xx,yy]=0.
	for xx=0,61 do for yy=0,61 do if sqrt(((xx-31.)^2.)  +  ((yy-31.)^2.)) gt annmode_inout[1] then imgsdi[xx,yy]=0.

writefits, folder+obj+'_sdiklip_conv'+suffix+'.fits', imgsdi, head

writefits, folder+obj+'_sdiklip_xyn_cube'+suffix+'.fits', sdis, head

path=folder
oldfiles=FILE_SEARCH(path,'*prelim.fits',COUNT=cc)
FILE_DELETE, oldfiles


endif

;ales_spec

;comb_mwc758

;oldnew=[ [[old]],[[new]] ]

;writefits,'/Volumes/Volumes/RAID36TB/LMIRCam/PDS201_4/reduced/PDS201_sdiklip_oldnew.fits',oldnew

end
