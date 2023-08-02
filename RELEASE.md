### Changes
- Removed the option to run at higher verbosities, will instead throw an error if maxima is far from atom.
- Added flag to pass to bca to control the distance at which the maxima distance error is thrown.
- Now runs with 1E-6 vacuum tolerance as default.
- Reduced memory usage of VoxelMap.
- Split out sum_bader_densities to calculate_bader_density and calculate_bader_volume_radius
