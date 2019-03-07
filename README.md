# immport
Manage Immport data for Hieracrchical Community Network

## *_converter.r 
### step 1
Pull datatype from each study from immport using ImmuneSpaceR
Data type: elisa, elispot, fcs, antibody (neut_ab_titer,hai)

### step 2 
Generate data matrix by taking median group by population and observation 

### step 3 
Generate annotation by subset on input file 

## organize by study.r 
Tidy up each dataset by study name 

## ispace_query.r 
Create matrix on study name and dataset , summarize on what data set are in each study

## biosample_processor 
Connect to sql database and pull biosample.txt for each study 
