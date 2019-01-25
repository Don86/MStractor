# MStractor
MStractor, an R workflow for non-targeted processing of LC-MS data

The MStractor workflow performs the following:
1.	Feature extraction (m/z-retention-time pair) using XCMS.
2.	Retention time alignment of the detected features across the samples composing the analytical set.
3.	Recognition and annotation of isotope clusters, fragments and charge states using CAMERA.
4.	Molecular feature filtering based on multiple criteria and conservative usage of intensity thresholds for maximum sensitivity.
5.	Creation of a data matrix summarizing the data processing results.

MStractor shows some additional features such as: 
1.	Parameterization based on user provided inputs obtained from instrument specifications and reference measurements.
2.	Graphical tools for real-time quality monitoring and optimization of the feature extraction process


