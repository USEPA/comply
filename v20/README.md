# RELEASE NOTES for COMPLY v2.0, 10/22/2026

More information about COMPLY is available at epa.gov/radiation/comply
***
New features, major updates, and bug fixes are described below compared to COMPLY v1.7.1. 

### New Features and Major Updates
* No changes to science or model inputs/parameters
* Code updated to modern Fortran convention (new version will be v2.0)
* User interface changed to a namelist-based (text user interface; TUI)
* Numerous binary input files converted to a single ascii format CSV file
* Standard output more comprehensive and added a new diagnostic output file
  + Standard output provides nuclide specific dose by exposure pathway
  + Diagnostic output provides nuclide specific air concentration and dose at evenly spaced receptors from 20 m to 50 km downwind
* Allow user input for annual rainfall rate
* Source code available along with scripts to compile on Windows, Mac, or Linux computing systems
  
### Bug Fixes
* Bug in wet deposition calculation which only impacts simulations using Level 4 with a wind rose.
  + Implications: Slightly higher dose predictions.


