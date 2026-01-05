FUNCTION GET_KLIP_BASIS, image_arr, k_klip
  compile_opt IDL2
  
  ; Function that returns the KL basis of input matrix
  ; IMAGE_ARR - Data matrix : N_obj (columns) X M_attrib (rows) ; i.e., column vectors
  ; K_KLIP - Integer number of KL Basis Vectors to keep
  
  ; Save initial k_klip state
  k_klipi = k_klip
  
  ; Get dimensions
  sz = size(image_arr)
  R_Dimen = sz[0]
  
  ; Efficient dimension handling
  if R_Dimen eq 1 then begin
    n_col = sz[1]
    n_row = 1
  endif else begin
    n_col = sz[1]
    n_row = sz[2]
  endelse
  
  ; Limit k_klip to number of available columns
  k_klip = k_klip < n_col
  
  ; Data Matrix "R"
  R = image_arr
  
  ; Constants
  TOLERANCE = 5.0E-20  ; Close enough to zero?
  
  ; Compute covariance matrix
  A = R # transpose(R)
  
  ; Get eigenvalues/vectors with EIGENQL
  eigenval = eigenql(A, EIGENVECTORS = eigenvect, /double)
  
  ; Process eigenvalues and normalize eigenvectors
  n_eig = n_elements(eigenval)
  for jj = 0, n_eig-1 do begin
    if abs(eigenval[jj]) lt TOLERANCE then eigenval[jj] = 0.0
    if eigenval[jj] gt 0 then begin
      eigenvect[*,jj] = eigenvect[*,jj]/sqrt(eigenval[jj]) ; normalize
    endif
  endfor
  
  ; Calculate new basis
  new_basis = matrix_multiply(eigenvect, R, /atranspose)
  
  ; Extract only the k_klip components we need
  out_sz = size(new_basis)
  new_basis_Dimen = out_sz[0]
  if new_basis_Dimen eq 1 then n_rows = 1 else n_rows = out_sz[2]
  
  kl_basis = fltarr(k_klip, n_rows)
  for jj = 0, k_klip-1 do begin 
    kl_basis[jj,*] = new_basis[jj,*]
  endfor
  
  ; Restore initial k_klip state
  k_klip = k_klipi
  
  print, 'created KLIP Basis.'
  return, kl_basis
END