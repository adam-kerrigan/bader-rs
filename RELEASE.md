## v0.5.0
### Features
- Added critical point handling.
- Atoms can now be bonded to images of themselves.
### Changes
- Voronoi shift has been changed to include the passing of periodic boundaries, a no check version has been added to perform as previously.
- The weight storage has changed adding EncodedWeight, EncodedAtom and EncodedImage.
- Added helper functions for parallelising.
- Removed crossbeam from dependencies.
- Removed rustc_hash from dependencies.
- Add own hash functions.
- Upped maximum distance flag to 1 Angstrom.
