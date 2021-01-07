## v0.3.0
### Changes
- Progress bars changed to not show by default, stops drawing of default bar
- Functions that create a progress bar (assign_atom, charge_sum) now take one as an argument
### Bug Fixes
- Overflow error when vacuum tolerance is so high that all charge is vacuum
## v0.2.3
### Changes
- Set up new target for releases that doesn't require GLibC
- Removed parking_lot::Mutex from the main program
- Set up Zenodo
### Bug Fixes
- Corrected the total of the Assigning to Atoms progress bar
## v0.2.2
### Bug Fixes
- Fixed SegFault at high thread count by pre-allocating weight_map ([issue: #19](https://github.com/kerrigoon/bader-rs/issues/19))
### Feature Changes
- Added a cap of 12 to the amount of threads distrubuted over by default
## v0.2.1
### Bug Fixes
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
