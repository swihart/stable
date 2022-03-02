# stable R package
Bruce Swihart  
Mar 2022

## Submission 1

  * refined the previous version's references to contemporary stable distribution 
  literature to help contextualize this package in the /man files and README


## Test environments
Local OS X: R version 4.1.2 (2021-11-01)
  * Platform: x86_64-apple-darwin17.0 (64-bit)
  * Running under: macOS Mojave 10.14.6
  
rhub::check(platform = "debian-clang-devel"): Debian Linux, R-devel, clang, ISO-8859-15 locale
rhub::check(platform = "windows-x86_64-devel"): Windows Server 2022, R-devel, 64 bit


## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.

# stable R package
Bruce Swihart  
Feb 2022

## Submission 1

  * http: -- > https: where appropriate
  * added references to contemporary stable distribution literature to help
  contextualize this package in the /man files and README
  * Removed imports: stabledist from DESCRIPTION
    

## Test environments
Local OS X: R version 4.1.2 (2021-11-01)
  * Platform: x86_64-apple-darwin17.0 (64-bit)
  * Running under: macOS Mojave 10.14.6
  
rhub::check(platform = "debian-clang-devel"): Debian Linux, R-devel, clang, ISO-8859-15 locale
rhub::check(platform = "windows-x86_64-devel"): Windows Server 2022, R-devel, 64 bit
rhub::check(platform = "ubuntu-rchk")

devtools::check_win_devel()
devtools::check_win_release()
https://win-builder.r-project.org/326a56Lp6aAx/

macOS M1 Builder: https://mac.r-project.org/macbuilder/results/1645028678-c805840b15690b3d/

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.


# stable R package
Bruce Swihart  
Feb 2019

## Re-Submission 1

  * fixed  All declared Imports should be used. NOTE.
  * updated stable.r with Roxygen comments

## Test environments
* local OS X install: R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
* Ubuntu 14.04.5 LTS (on travis-ci): R version 3.5.2 (2017-01-27)
* win-builder: R Under development (unstable) (2019-01-31 r76038)

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.

# stable R package
Bruce Swihart  
May 2018

## Submission 1 of v1.1.3:
  Fixed note; added functions.

## Test environments
* local OS X install: R version 3.4.1 (2017-06-30)
* ubuntu 12.04 (on travis-ci): R version 3.5.0 (2017-01-27)
* win-builder: R Under development (unstable) (2018-05-20 r74747) -- "Unsuffered Consequences"

## R CMD check results
There were no ERRORs or WARNINGs or Notes. 

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Miscellaneous
Added some functions.