## Bader-rs v0.2.0
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
