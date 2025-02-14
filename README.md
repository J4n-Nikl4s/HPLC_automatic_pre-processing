Research data management (RDM) is important to create reproducible and reliable data but can be a time-consuming. In this project, we established an RDM pipeline (Figure 1 A) for HPLC data created by the Agilent HPLC-software (OpenLab CDS 2.7) (Figure 1 B), using the Python package [chromatopy](https://github.com/FAIRChemistry/chromatopy) developed by Max Häußler. This helps to speed up and standardize RDM for the HPLC analysis. The data are saved in a machine-readable XML format called [EnzymeML](https://enzymeml.github.io/tools/).  We created this project to help us with the analysis of enzyme and enzyme cascade reactions but could be used in theory for other use cases.

![image](https://github.com/user-attachments/assets/14bd9ae4-0af8-4848-8ec3-3c960c839fcf)


Figure 1 Schematic overview of the research data management piepline 

## How to use
For information how to install visite the [chromatopy documentation]([https://github.com/FAIRChemistry/chromatopy](https://fairchemistry.github.io/chromatopy/#installation)) page
1. Data structure
   
   For the input of the data structure, shut ordered like represented in Figure 2.
   
   The calibration **folder** shut include the molecule name and the [PubChem](https://pubchem.ncbi.nlm.nih.gov/) Compound CID. The calibration data **files** shut include the concentration of the compound and the unit.
   The time sample data **files** shut include the time the sample was taken and the unit.
   
   ![image](https://github.com/user-attachments/assets/3340c58b-4406-4947-9016-6440842a9170)

   
   Figure 2 Representation of the data structure

3. From the calibration data, the script HPLC_creation_molecuels creates a virtual molecule.
4. Depending on how much data has to be analyzed, two different scripts were created.

   HPLC_single_batch
   In this case, there are only a few experiments with always the same parameter or the same enzyme.

   HPLC_multi_batch
   In this case, there are many experiments with different parameters, such as temperature, pH, substrate combinations, etc. 


3. Both scripts are creating an EnzymeML_file.json and  concentration_time_courses.png. 
