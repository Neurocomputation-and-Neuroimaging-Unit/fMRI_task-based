# fMRI_task-based

This a pipeline dedicated to the conversion, preprocessing, univariate and multivariate analysis 
on subject and group level.
It includes 5 steps:

  - A: convert dicom data to 4d-nifti data in bids structure
    adjust for nature of your data in A0
  - B: preprocess data 
    possible preprocessing steps (choose in B0) are:
      - Segemtation
      - Slice-time correction
      - Realignment
      - Coregistration (estimate, optional: reslice)
      - Normalization
      - Scrubbing
      - Smoothing
      - Detrending
  - C: univariate analysis
    adjust for your design and select data in C0
      - reading onsets (bids-conform tsv per run on subject-level required!)
      - glm specification
      - contrasts specification
      - 2nd level: 1 folder per contrast
      - 2nd level: flex fact
  - D: multivariate analysis
    adjust for your design in D0, possible steps include:
      - reading onsets (bids-conform tsv per run on subject-level required!)
      - fir specification
      - glm specification
      - regression (SVR) or classification (SVM), searchlight or ROI-based
      - averaging for convenience
      - Normalization (eg. for accuracy maps or zscore maps)
      - Smoothing
  - E: flexible factorial design for decoding results
  
  Ideas for updates are welcome!
