# SUNTANS

The Stanford unstructured-grid, nonhydrostatic, parallel coastal ocean model. For simulation of nonhydrostatic flows at high resolution in estuaries and coastal seas. Requires a grid generator and ParMETIS (if run in parallel).

Created by Oliver Fringer

<img src="./guide/images/galviston.png" height="350" title="Unstructured grid of Galveston Bay, TX (Matt Rayson)"> <img src="./guide/images/snoho3.jpg" height="350" title="Near Surface Salinity in the Snohomish River (Bing Wang)">

### Quickstart Guide
Download the main source and examples from the GitHub repository here:
https://github.com/ofringer/suntans/archive/master.zip

Follow these steps:

1.  `unzip suntans-master.zip`
2.  `cd suntans-master/examples/cylinder`
3.  `make test`
4.  `cd ../../main`
5.  `make sunplot`
6.  `./sunplot --datadir=../examples/cylinder/data`
