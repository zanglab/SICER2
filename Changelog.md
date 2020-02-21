# Changelog
All notable changes to this project will be documented in this file.


## [1.0.2] - 2020-02-21
### Added
- Changelog file to keep track of changes to this project.

### Fixed
- Fixed bug with outputting to an user-specified output directory.
- Fixed bug in `find_islands_in_pr` so that empty island lists do not trigger `IndexError`. 

## [1.0.1] - 2019-12-28
### Fixed
- Fixed issue with calling BEDTools using `subprocess.call`. Instead, use `subprocess.Popen` to enable piping.
- Fixed method for getting CPU count
    
## [1.0.0] - 2019-06-14
### Added
- New version of SICER for improved user-friendliness and parallelization support
- Multiprocessing based on chromosomes
- RECOGNICER algorithm added
- Uses BEDTools to convert BAM file to BED automatically