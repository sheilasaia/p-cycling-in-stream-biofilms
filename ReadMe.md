# "paper-p-cycling-in-stream-biofilms" ReadMe #

This ReadMe.md file was generated on 20170112 by Sheila Saia.

This GitHub repository was created to provide access to collected data, analysis code, and other information associated with the paper by Saia et al. titled 'Evidence for polyphosphate accumulating organism (PAO)-mediated phosphorus cycling in stream biofilms under alternating aerobic/anaerobic conditions' in Freshwater Science.

## General Information ##

**Title of Dataset**<br>
"paper-p-cycling-in-stream-biofilms"

**Contact Information**<br>
Name: Sheila Saia<br>
Institution: Cornell University<br>
Address: B62 Riley-Robb Hall, Ithaca, NY 14853<br>
Email: sms493@cornell.edu<br>

**Date of data collection**<br>
These data were collected during a laboratory experiment from 20151003 to 20151005.

**Geographic location of data collection**<br>
These data were collected in Cascadilla Creek, Ithaca, NY as well as in the Soil & Water Lab at Cornell University in Ithaca, NY.

**Information about funding sources that supported the collection of the data**<br>
Sheila Saia was supported by a Cornell University CALS Land Grant Fellowship and USEPA STAR Fellowship. This project was supported by a USDA NIFA grant #2014-67019-21636.

## Sharing & Access Information ##

**Licenses/restrictions placed on the data**<br>
Please use and distribute according to CC-BY v4.0. For a human readible version of this license visit [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/) .

**Links to publications that cite or use the data**<br>
As of 20161208 there are no other publications that cite or use these data.

**Links to other publicly accessible locations of the data**<br>
This dataset and associated R code are available at https://github.com/sheilasaia/paper-p-cycling-in-stream-biofilms (DOI:XXXXX). The associated publication is available via Freshwater Science (DOI:XXXXX).

**Links/relationships to ancillary data sets**<br>
There are no links to or relationships with other ancillary data sets.

**Data derived from another source**<br>
Data was not derived from another source.

**Recommended citation for the data**<br>
Saia, S. M., P. J. Sullivan, J. M. Regan, H. J. Carrick, A. R. Buda, N. A. Locke, M. T. Walter. 2017. Evidence for polyphosphate accumulating organism (PAO)-mediated phosphorus cycling in stream biofilms under alternating aerobic/anaerobic conditions. Freshwater Science. XXXXX(XXXXX):XXXXX-XXXXX.

## Data & File Overview ##

**File List**<br>
Filename: allTUBdata\_oct2014\_forPaper.txt <br>
Short description: This text file includes water quality related data for this experiment including: phosphate concentrations, FeII concentrations, and cation concentrations. It also includes environmental variables measured during the experiment including pH, dissolved oxygen, and temperature.<br>

Filename: cellCounts\_oct2014\_forPaper.txt <br>
Short description: This text file includes the cell counts (DAPI-DNA and DAPI-polyphosphate) taken from both treatments at the end of the experiment.<br>

Filename: PPextAll\_oct2014\_forPaper.txt<br>
Short description: This text file includes the results of the polyphosphate biofilm extractions at the start and ends of the experiment.<br>

Filename: TPextAll\_oct2014\_forPaper.txt<br>
Short description: This text file includes the results of the total phosphate biofilm extractions at the start and ends of the experiment.<br>

Filename: oct2014experiment\_script\_forPaper\_final.R<br>
Short description: This R script includes code for all data analysis (apart from calibration of phosphate, FeII, and ICP-MS results) including the statistical analysis and data visualization used in the Freshwater Science journal article associated with these data.<br>

Filename: oct2014experiment\_script\_forPaper\_final.Rmd<br>
Short description: This file is the same as the oct2014experiment\_script\_forPaper\_final.R file but in the R markdown format.<br>

**Relationship Between Files**<br>
The text files listed above (allTUBdata\_oct2014\_forPaper.txt, cellCounts\_oct2014\_forPaper.txt, PPextAll\_oct2014\_forPaper.txt, TPextAll\_oct2014\_forPaper.txt) are all required for running the R script called oct2014experiment\_script\_forPaper\_final.R. The oct2014experiment\_script\_forPaper\_final.Rmd file is identical to oct2014experiment\_script\_forPaper\_final.R, but in the R mark-up format.

**Raw Data**<br>
This repository also contains the raw data that was used compiled into the files listed above. Raw data can be found in the directory called raw\_data. Sub-directories within this main directory include:<br>

Directory name: cell\_counts <br>
Short description: This directory contains microscope images (DAPI-DNA and DAPI-polyphosphate (polyP)) for all treatment replicate views as well as the grid drawing used for manual cell counts. Methods are described in the associated Freshwater Science journal article.<br>

Directory name: feII\_analysis <br>
Short description: This directory contains the raw and processed FeII data obtained using the Ferrozine assay and Tecan plate reader described in the methods section of the associated Freshwater Science journal article.<br>

Directory name: icp\_analysis <br>
Short description: This directory includes the raw ICP-MS data from this experiment and described in the methods section of the associated Freshwater Science journal article.<br>

Directory name: srp\_analysis <br>
Short description: This folder contains the raw phosphate (aka srp) data from this experiment (rawPdata\_20141013.xlsx), the R script used to calibrate these data (P\_calibration\_code\_20141013.R), text files that are used as inputs for the R script (pCalib\_, pCheck\_, pData\_, pDI\_), text file outputs of the R script (pNewData\_), and a file that combines these inputs and outputs (Pdata\_20141013.xlsx).

Filename: oct2014experiment\_rawData\_forPaper.xlsx <br>
Short description: This file has all the processed data in one place and is the source of the text files mentioned in the 'File List' section of this ReadMe file (see tabs with the same names). Briefly, the 'schedule' tab defines the sampling timing and protocol for this experiment as well as the treatment labels. These labels are also explained in the associated Freshwater Science journal article.<br>

**Additional related data collected that was not included in the current data package:**<br>
All data is included in this package.

**Are there multiple versions of the dataset?**<br>
No, there are no other versions of this dataset.

## Methodological Information ##

**Description of methods used for collection/generation of data:**<br>
See the associated Freshwater Science journal article for a full description of the methods used to collect and analyze these data.

**Methods for processing the data:**<br>
See the R scripts in this repository as well as the associated Freshwater Science journal article for a full description of the methods used to collect and analyze these data.

**Instrument- or software-specific information needed to interpret the data:**<br>
The latest version of the R programming lanugage is required to run the R scripts in this repository. R can be downloaded  for free here: [https://www.r-project.org/](https://www.r-project.org/). Microsoft Excel is required to open .xlsx files.

**Standards and calibration information, if appropriate:**<br>
Information on calibrations are included in the 'Raw Data' section of this ReadMe file.

**Environmental/experimental conditions:**<br>
See the associated Freshwater Science journal article for a full description of the environmental and experimental conditions used while collecting these data.

**Describe any quality-assurance procedures performed on the data:**<br>
We double checked manually entered data and plotted data in R (scripts not encluded) to ensure that errors were not made when manually entering data.

**People involved with sample collection, processing, analysis and/or submission:**<br>
See the associated Freshwater Science journal article for a full description of author contributions and acknowledgments.

## Data-Specific Information For: allTUBdata\_oct2014\_forPaper.txt ##

**Variable list**<br>
'HoursFromStart' - Hours from start of experiment.<br>
'SampleNum' - Sample number.<br>
'TubID' - Tub treatment identifier where T1 represents alternating anaerobic/aerobic treatment and T2 represents control treatment that was always aerobic.<br>
'Condition' - Identifier to explain specific condition of treatment for a given sample number where 'air' means the treatment was being bubbled with air and where 'n2' means the treatment was being bubbled with a mixed anaerobic gas (80% N2:20% CO2 gas).<br>
'pH' - pH of overlying water for a given sample.<br>
'DOmgL' - Dissolved oxygen concentration (mg/L) of overlying water for a given sample.<br>
'AvgTempC' - Average temperature concentration (degrees C) of the overlying water for a given sample.<br>
'AvgSRPppm' - Average soluble reactive phosphorus (i.e. phosphate) concentration (ppm) of the overlying water for a given sample as analyzed by the molybdenum blue method.<br>
'StdSRPppm' - Standard deviation of the average soluble reactive phosphorus (i.e. phosphate as P) concentration (ppm) of the overlying water for a given sample as analyzed by the molybdenum blue method.<br>
'TotalPmgL' - Total phosphorus concentration (mg/L of the overlying water for a given sample as analyzed by ICP-MS.<br>
'Fe2ppm' - Iron II concentration (ppm) of the overlying water for a given sample as analyzed by ferrozine assay.<br>
'TotalFemgL' - Total iron concentration (mg/L) of the overlying water for a given sample as analyzed by ICP-MS.<br>
'TotalCamgL' - Total calcium concentration (mg/L) of the overlying water for a given sample as analyzed by ICP-MS.<br>
'TotalSmgL' - Total sulfer concentration (mg/L) of the overlying water for a given sample as analyzed by ICP-MS.<br>
'TotalKmgL' - Total potassium concentration (mg/L) of the overlying water for a given sample as analyzed by ICP-MS.<br>
'TotalMgmgL' - Total magnesium concentration (mg/L) of the overlying water for a given sample as analyzed by ICP-MS.<br>
'TotalMnmgL' - Total mangenese concentration (mg/L) of the overlying water for a given sample as analyzed by ICP-MS.<br>

**Missing data codes**<br>
NA - Data was below detection limit of machine.

## Data-Specific Information For: cellCounts\_oct2014\_forPaper.txt ##

**Variable list**<br>
'SampleID' - Indentfier to explain the specific treatment at the end of the experiment for a given field of view where T1 represents anaerobic treatment and T2 represents control treatment that was aerobic.<br>
'SampleRep' - Number of replicate for a particluar treatment.<br>
'ViewNumber' - Number of field of view for a particluar sample.<br>
'DNACounts' - Number of cells fluorescing under DAPI-DNA filter set for a particular field of view.<br>
'PolyPCounts' - Number of cells fluorescing under DAPI-polyphosphate filter set for a particular field of view.<br>
'Dillution' - Dilution strength of a given sample.<br>
'DNAcellsmL' - Total number of cells per mL of sample.<br>
'PolyPcellsmL' - Total number of cells with stored polyphosphate granules per mL of sample.<br>
'PerPolyP' - Percent number of cells with stored polyphosphate granules for a given sample.<br>

**Missing data codes**<br>
NA - Cells in view had no polyP granules.

## Data-Specific Information For: PPextAll\_oct2014\_forPaper.txt ##

**Variable list**<br>
'Num' - Row number.<br>
'SampleID' - Indentfier to explain the specific treatment at the end of the experiment for a given field of view where T1 represents anaerobic treatment and T2 represents control treatment that was aerobic.<br>
'Replicate' - Replicate number.<br>
'Extraction' - Type of extraction. PP for polyphosphate or TP for total phosphate.<br>
'Period' - Indication whether biofilm sample was taken for analysis at the start or end of the experiment.
'myPppmFix' - Calibrated polyphosphate concentration (ppm) in the digested biofilm (as soluble reactive phosphorus or phosphate as P). Analyzed by the molybdenum blue method.<br>
'wetBFg' - Mass of wet biofilm added to the digestion in grams.<br>
'VolAddedmL' - Volume of solution added to the digestion in mL.<br>
'AvgSAm2' - Average surface area of the cobbles for that particular treatment and period in m^2.<br>

**Missing data codes**<br>
No missing data codes.

## Data-Specific Information For: TPextAll\_oct2014\_forPaper.txt ##

**Variable list**<br>
'Num' - Row number.<br>
'SampleID' - Indentfier to explain the specific treatment at the end of the experiment for a given field of view where T1 represents anaerobic treatment and T2 represents control treatment that was aerobic.<br>
'Replicate' - Replicate number.<br>
'Extraction' - Type of extraction. PP for polyphosphate or TP for total phosphate.<br>
'Period' - Indication whether biofilm sample was taken for analysis at the start or end of the experiment.
'myPppmFix' - Calibrated total phosphorus concentration (ppm) in the digested biofilm (as soluble reactive phosphorus or phosphate as P). Analyzed by the molybdenum blue method.<br>
'wetBFg' - Mass of wet biofilm added to the digestion in grams.<br>
'VolAddedmL' - Volume of solution added to the digestion in mL.<br>
'AvgSAm2' - Average surface area of the cobbles for that particular treatment and period in m^2.<br>

**Missing data codes**<br>
No missing data codes.