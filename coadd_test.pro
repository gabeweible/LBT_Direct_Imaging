pro coadd_test
  coadds = [8, 10, 20, 25, 40, 50, 100, 200]; Allowed coadds that actually can run on my machine
  foreach coadd, coadds do hii1348_pipeline, neg_inj=0, pre_inj_stuff=1, coadd=coadd
end
