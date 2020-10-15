## v0.3.0
### Bug fixes
- Added a lock to maxima_get() in VoxelMap and made a maxima_non_blocking_get(), unsure if this would ever be a problem due to the lock on index.pop() but better safe than sorry.
### Documentation Changes
- More of the crate documented
- Documentation tests added for all partitioning methods and for using weight_store in VoxelMap
## v0.2.0
### New Features
- Added spin flag for allowing cubes have spin and density output
- Complete revamp weight method, now very fast and scales well
### Removed Features
- No longer able to apply weighting of boundaries to all methods, just the weight method
### UI Changes
- '--weight, -w' now controls the weight tolerance allowing extremely small contributions to be discarded
- '--spin, -s' has been added for spin output on cube files
### Library Changes
- I/O now has a trait for standardising implementation of new file types
- I/O modules now contribute to the formating of ACF and BCF files (maybe units in future)
- VoxelMap now controls the population and processing of the voxel maps
- Custom Lock introduced to speed up threading
## v0.1.1
### Bug Fixes
- Volume weighting logic error ([issue: #5](https://github.com/kerrigoon/bader-rs/issues/5))
### Cosmetic Changes
- Standardised the method description in help information
- Changed the BCF.dat file to include atom number
