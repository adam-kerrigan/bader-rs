### Changes
- Removed the need of passing a density to the calculate_bader_volume_radius function.
- Changed the name of calculate_bader_volume_radius to calculate_bader_volumes_and_radii.
- Changed AtomVoxelMap to VoxelMap as there are no longer two VoxelMap variants.
- removed the VoxelMap triat.
- Changed the name of VoxelMap.boundary_iter() to VoxelMap.weight_iter().
- Added VoxelMap.maxima_len() and VoxelMap.weight_len().
- Removed anyhow for the error management.
- Changed the return of invert_lattice to Option as there is only one way it can fail.
