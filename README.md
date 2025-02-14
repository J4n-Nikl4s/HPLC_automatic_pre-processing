Research data management (RDM) is important to create reproducible and reliable data and can be a time-consuming. In this project, we established an RDM pipeline for HPLC data created by the Agilent HPLC-software (OpenLab CDS 2.7), using the Python package [chromatopy](https://github.com/FAIRChemistry/chromatopy) developed by Max Häußler. This helps to speed up and standardize RDM for the HPLC analysis. The data are saved in a machine-readable XML format called [EnzymeML](https://enzymeml.github.io/tools/).  We created this project to help us with the analysis of enzyme and enzyme cascade reactions but could be used in theory for other use cases.

1. From the calibration data, the script HPLC_creation_molecuels creates a virtual molecule.
2. Depending on how much data has to be analyzed, two different scripts were created. 

### HPLC_singel_batch
In this case, there are only a few experiments with always the same parameter or the same enzyme.

### HPLC_multi_batch
In this case, there are many experiments with different parameters, such as temperature, pH, substrate combinations, etc. 
3. Both scripts are creating an EnzymeML_file.json and a concentration_time_courses.png. 
