# stable R package
Bruce Swihart  
December 2016

## Test environments
* local OS X install: R version 3.3.2 (2016-10-31)
* ubuntu 12.04 (on travis-ci): R version 3.3.1 (2016-06-21)
* win-builder: R Under development (unstable) (2016-12-11 r71774)

## R CMD check results
There were no ERRORs or WARNINGs or Notes. 

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Miscellaneous
Prof Ripley emailed me and CRAN and asked me to fix a Valgrind error, which this patch/resubmission does.

Email:
On Wed, Dec 14, 2016 at 4:24 AM, Prof Brian Ripley <ripley@stats.ox.ac.uk> wrote:

    On 11/12/2016 10:16, Prof Brian Ripley wrote:

        See https://cran.r-project.org/web/checks/check_results_stable.html .

        This clearly attempting to access the non-existent element 5 of tab1[5].

        The background to detecting the fatal SAN errors is in 'Writing R
        Extensions': please correct ASAP.

    Unfortunately you have replaced this by use of uninitialized memory: see the 'valgrind' link now on that page.

