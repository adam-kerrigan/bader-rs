## v0.5.0
### Features
- Added critical point handling.
- Atoms can now be bonded to images of themselves.
### Changes
- Voronoi shift has been changed to include the passing of periodic boundaries, a no check version has been added to perform as previously.
## v0.4.8
### Changes
- Opt -> Arg in arguments.
- added a parse_filetype function.
## v0.4.7
### Small things
- Took VoxelMap out of the box in bca.
## v0.4.6
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
## v0.4.5
### Changes
- Removed the need of passing a density to the calculate_bader_volume_radius function.
- Changed the name of calculate_bader_volume_radius to calculate_bader_volumes_and_radii.
- Changed AtomVoxelMap to VoxelMap as there are no longer two VoxelMap variants.
- removed the VoxelMap triat.
- Changed the name of VoxelMap.boundary_iter() to VoxelMap.weight_iter().
- Added VoxelMap.maxima_len() and VoxelMap.weight_len().
- Removed anyhow for the error management.
- Changed the return of invert_lattice to Option as there is only one way it can fail.
## v0.4.4
### Bug fix
- Removed index deletion in maxima finding.
## v0.4.3
### Changes
- Removed the option to run at higher verbosities, will instead throw an error if maxima is far from atom.
- Added flag to pass to bca to control the distance at which the maxima distance error is thrown.
- Now runs with 1E-6 vacuum tolerance as default.
- Reduced memory usage of VoxelMap.
- Split out sum_bader_densities to calculate_bader_density and calculate_bader_volume_radius
## v0.4.2
### Changes
- Updated dependancies and fixed the breaking changes associated with them.
- Bumped the minimum rust version required.
## v0.4.1
### Changes
- Moving to allow python bindings by separating functions and making structs more streamlined.
- Switch the entire analysis section to functions rather than a struct.
- Threaded charge summing, assigning maxima to atoms and the new maxima finding function.
- Lots of moving around of functions and changing outcomes, i.e, to_cartesian now returns cartesian coordinates.
- Updated to Clap v3.
- Created a VoxelMap Trait.
### Features
- Added a nearest neighbour function.
## v0.4.0
### Changes
- VoxelMap now handles the running of the bader calculation, using VoxelMap::calc().
## v0.3.2
### Features
- Writing of the charge density is now suppported.
### Changes
- Changed how the maxima and weights are stored for the boundary voxels.
- Memory optimisations ([issue: #30](https://github.com/adam-kerrigan/bader-rs/issues/30)).
## v0.3.1
### Changes
- Progress bars changed to not show by default, stops drawing of default bar.
- Functions that create a progress bar (assign_atom, charge_sum) now take one as an argument.
- Moved analysis functions (assign_atom, charge_sum, atoms_charge_sum) to own structure.
- Simplified VoxelMap to remove the weight index map. Weight indices are stored as negative numbers.
- Removed Rayon ([issue: #25](https://github.com/kerrigoon/bader-rs/issues/25)).
- Moved writing output to own module in anticipation of removing prettytable-rs.
- Dropped ongrid and neargrid, weight method fast and superior.
- github username change so updated all the links.
- Density now no longer contains the density and so has been renamed Grid.
- Crossbeam is used for threading scopes.
- Changed program name from bader to bca (Bader charge analysis).
- Removed Prettytables.
### Bug Fixes
- Overflow error when vacuum tolerance is so high that all charge is vacuum.
## v0.3.0 - (Yanked)
## v0.2.3
### Changes
- Set up new target for releases that doesn't require GLibC.
- Removed parking_lot::Mutex from the main program.
- Set up Zenodo
### Bug Fixes
- Corrected the total of the Assigning to Atoms progress bar.
## v0.2.2
### Bug Fixes
- Fixed SegFault at high thread count by pre-allocating weight_map ([issue: #19](https://github.com/kerrigoon/bader-rs/issues/19)).
### Feature Changes
- Added a cap of 12 to the amount of threads distrubuted over by default.
## v0.2.1
### Bug Fixes
- Added a lock to maxima_get() in VoxelMap and made a maxima_non_blocking_get(), unsure if this would ever be a problem due to the lock on index.pop() but better safe than sorry.
### Documentation Changes
- More of the crate documented
- Documentation tests added for all partitioning methods and for using weight_store in VoxelMap.
## v0.2.0
### New Features
- Added spin flag for allowing cubes have spin and density output.
- Complete revamp weight method, now very fast and scales well.
### Removed Features
- No longer able to apply weighting of boundaries to all methods, just the weight method.
### UI Changes
- '--weight, -w' now controls the weight tolerance allowing extremely small contributions to be discarded.
- '--spin, -s' has been added for spin output on cube files.
### Library Changes
- I/O now has a trait for standardising implementation of new file types.
- I/O modules now contribute to the formating of ACF and BCF files (maybe units in future).
- VoxelMap now controls the population and processing of the voxel maps.
- Custom Lock introduced to speed up threading.
## v0.1.1
### Bug Fixes
- Volume weighting logic error ([issue: #5](https://github.com/kerrigoon/bader-rs/issues/5)).
### Cosmetic Changes
- Standardised the method description in help information.
- Changed the BCF.dat file to include atom number.
