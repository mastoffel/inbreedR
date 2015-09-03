## Test environments
* local Windows 7, R 3.2.2
* ubuntu (on travis-ci), R 3.2.2
* win-builder (release)

## R CMD check results
There were no ERRORs or WARNINGs 

There was 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Martin A. Stoffel <martin.adam.stoffel@gmail.com>'
New submission


## Downstream dependencies
I have also run R CMD check on downstream dependencies of inbreedR.
All packages that I could install passed.
    
## Resubmission
This is a resubmission. In this version I have:

* Changed the Description field to not start with the package name.

* Made the title shorter and removed acronyms

* Assigned several functions to the stats or graphics package.