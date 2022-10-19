pro create_test
compile_opt idl2

raw_fits_path = '~/OneDrive/Research/HII1348/HII1348/raw'
print, 'Searching for FITS files in', raw_fits_path, '...'
files = FILE_SEARCH(raw_fits_path, '*.fits', COUNT=filecount)
print, 'Found ', filecount, ' FITS files!'
;
;obj_error = 0
;lmir_error = 0
;dit_error = 0
size_error = 0

;wrong_names = list()
;wrong_lmir = list()
;wrong_dit = list()
wrong_size = list()

start_frame = 0
for ii = start_frame, filecount-1 do begin
   print, 'File index', ii, '/', filecount-1
   frame = readfits(files[ii], head)
;   if strcompress(fxpar(head, 'OBJNAME'), /rem) ne 'HII1348' then begin
;      wrong_names.Add, ii
;      obj_error = 1
;   endif
;   if fxpar(head, 'LMIR_FW2') eq 'ND2.0-T1' then begin
;      wrong_lmir.Add, ii
;      lmir_error = 1
;   endif
;   dit = fxpar(head,'ITIME')
;   if not dit/906.51 gt 0.95 and dit/906.51 lt 1.05 then begin
;      wrong_dit.Add, ii
;      dit_error = 1
;   endif
   frame = reform(frame[*,*,1] - frame[*,*,0])
   if (size(frame))[2] gt 2000 then begin
      wrong_size.Add, ii
      size_error = 1
   endif
endfor

;if obj_error then print, 'You have an object name problem...'
;help, wrong_names
;help, obj_error
;
;if lmir_error then print, 'You have a LMIR_FW2 problem...'
;help, wrong_lmir
;help, lmir_error
;
;if dit_error then print, 'You have a dit problem...'
;help, wrong_dit
;help, dit_error

if size_error then print, 'You have a size problem...'
print, wrong_size
print, size_error

end