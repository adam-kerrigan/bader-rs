### Changes
- Progress bars changed to not show by default, stops drawing of default bar
- Functions that create a progress bar (assign_atom, charge_sum) now take one as an argument
- Moved analysis functions (assign_atom, charge_sum, atoms_charge_sum) to own structure
- Simplified VoxelMap to remove the weight index map. Weight indices are stored as negative numbers
- Removed Rayon ([issue: #25](https://github.com/kerrigoon/bader-rs/issues/25))
- Moved writing output to own module in anticipation of removing prettytable-rs
- Dropped ongrid and neargrid, weight method fast and superior
- github username change so updated all the links
- Density now no longer contains the density and so has been renamed Grid
- Crossbeam is used for threading scopes
- Changed program name from bader to bca (Bader charge analysis)
- Removed Prettytables
### Bug Fixes
- Overflow error when vacuum tolerance is so high that all charge is vacuum
