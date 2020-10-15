## Bader-rs v0.2.1
### Bug fixes
- Added a lock to maxima_get() in VoxelMap and made a maxima_non_blocking_get(), unsure if this would ever be a problem due to the lock on index.pop() but better safe than sorry.
### Documentation Changes
- More of the crate documented
- Documentation tests added for all partitioning methods and for using weight_store in VoxelMap
