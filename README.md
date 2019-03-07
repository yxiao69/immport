# immport
Manage immport data for Hieracrchical Community Network

## *_converter.r 
### step 1
pull datatype from each study from immport using ImmuneSpaceR
data type: elisa, elispot, fcs, antibody (neut_ab_titer,hai)

### step 2 
generate data matrix by taking median group by population and observation 

### step 3 
generate annotation by subset on input file 

## organize by study.r 
tidy up each dataset by study name 

## ispace_query.r 
create matrix on study name and dataset , summarize on what data set are in each study

## biosample_processor 
connect to sql database and pull biosample.txt for each study 
