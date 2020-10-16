## Bader-rs v0.2.2
### Bug Fixes
- Fixed SegFault at high thread count by pre-allocating weight_map ([issue: #19](https://github.com/kerrigoon/bader-rs/issues/19))
### Feature Changes
- Added a cap of 12 to the amount of threads distrubuted over by default
