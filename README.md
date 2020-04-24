# Grass Phenology

<h2> DATA </h2>
<h3> <a href="https://github.com/bnetobradley/grassphenology">This github repository</a> contains the raw data and code used to explore the phylogenetic structure of phenology in grass species.</h3>

<h4> In the /data folder of this repository, there are raw data files, these data correspond to specimens housed at the UBC Herbarium, and climate data derived from these. All specimens were scored for phenology (here defined as the presence of anthers), using the aid of a disecting microscope, between September 2017 and March 2018. Climate data was downloaded in March 2018. </h4>

<h3>1. Herbarium data</h3>
<h5> Files: Find these under /data/ </h5>
<h5> Data Specifics: This copy of the UBC herbarium database was downloaded on the 13th of June 2018. </h5>

<h3>2. Climate data</h3>
<h5> Files: Find year by year records in /data/files_from_climate_na </h5>
<h5> Data Specifics: These data contain two ID columns, theh first corresponds to the herbarium specimen accession number from which this climate data is generated, the second is an automatically generated number that is essentially void of useful information. </h5>

<h2> CODE </h2>
<h4> In the /R folder of this repository, there are multiple files involved in cleaning and manipulating the datasets used. Below is a roadmap. </h4>

<h3>1. Cleaning: </h3>

<h3>2. Merging datasets: </h3> 
<h4>In this script we merge herbarium records with climate data derived from each record's latitute and longitude. We also clean up binomials, dates and remove species with insufficient flowering records (ie. fewer than 10.) </h4>

<h3>3. Variable selection: </h3>

<h3>4. Phylogenetic analyses: </h3>

