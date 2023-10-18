### Features
- Added a method for calculating the Laplacian at a point.
- Added way to calculate the error in the partitioning from the Laplacian.
### Changes
- Voronoi now stores the volume of the Voronoi cell.
- Removed clap as a dependancy.
- Changed the flags for file type to -f --file_type from -t --type.
- Changed the short flag for threads to -t from -J.
- Removed indicatif and atomic-counter as dependancies.
- Removed regex from dependancies.
- Progress bars are now created inside functions and whether they are shown is optional.
- Added a new silent flag: -x --silent.
