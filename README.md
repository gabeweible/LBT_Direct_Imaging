# LBT_Direct_Imaging_Pipeline
Various IDL procedures as part of a direct imaging data reduction and analysis pipeline. This pipeline is based off of work by Dr. Kevin Wagner and Dr. Daniel Apai of Steward Observatory, but has been heavily modified, refactored, and extended to fit my work. This is a work-in-progress and not intended to be a standalone package for other work, at the moment. Currently, I am using this to reduce LBTI direct imaging data of HII 1348. HII1348_pipeline.pro is the main pipeline procedure, and calls the other procedures in the pipeline directory as steps in the reduction (coadding, creating a data cube, bad pixel removal, sky background subtraction, KLIP and/or ADI, etc.)

![LBT-Enhanced](https://user-images.githubusercontent.com/116225423/217911189-497e39b6-cfd4-409a-b53e-ab145bfdfece.jpg)
Image from NASA/JPL-Caltech
