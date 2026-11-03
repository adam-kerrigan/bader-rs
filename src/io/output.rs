use crate::analysis::Partition;
use crate::atoms::Atoms;
use crate::critical::{CriticalPoint, CriticalPoints};
use crate::grid::Grid;
use crate::io::{FileFormat, FileType};
use crate::methods::laplacian;
use crate::voxel_map::EncodedAtom;
use core::fmt;
use std::fs::File;
use std::io::Write;

/// The names of the densities in a calculation: the charge followed by the
/// spin components.
pub fn density_names(densities: usize) -> &'static [&'static str] {
    match densities {
        0 | 1 => &["Charge"],
        2 => &["Charge", "Spin"],
        _ => &["Charge", "Spin X", "Spin Y", "Spin Z"],
    }
}

/// A row of the critical point table of an atom.
struct CriticalPointRow {
    kind: &'static str,
    position: [String; 3],
    density: String,
    laplacian: String,
    members: String,
}

/// The critical point information of every atom, formatted as markdown
/// tables by its [`Display`](fmt::Display) implementation.
pub struct CriticalPointOutput {
    summary: String,
    positions: Vec<[f64; 3]>,
    coordination: Vec<usize>,
    rows: Vec<Vec<CriticalPointRow>>,
}

impl CriticalPointOutput {
    /// Pivots the critical points onto the atoms they involve.
    ///
    /// # Arguments
    /// * `atoms`: The atoms of the system.
    /// * `critical_points`: The pruned critical points.
    /// * `density`: The reference density.
    /// * `grid`: The grid the density is stored on.
    /// * `file_type`: Used to format the coordinates for the file format.
    pub fn new(
        atoms: &Atoms,
        critical_points: &CriticalPoints,
        density: &[f64],
        grid: &Grid,
        file_type: FileType,
    ) -> Self {
        let atom_len = atoms.positions.len();
        let mut rows: Vec<Vec<CriticalPointRow>> =
            (0..atom_len).map(|_| Vec::new()).collect();
        let mut coordination = vec![0; atom_len];
        for (kind, cps) in [
            ("Nucleus", &critical_points.nuclei),
            ("Bond", &critical_points.bonds),
            ("Ring", &critical_points.rings),
            ("Cage", &critical_points.cages),
        ] {
            // sort so the output does not depend on thread scheduling
            let mut cps = cps.iter().collect::<Vec<_>>();
            cps.sort_unstable_by_key(|cp| cp.position);
            for cp in cps {
                let position = file_type
                    .coordinate_format(grid.to_cartesian(cp.position))
                    .map(|p| format!("{:.6}", p));
                let density_text =
                    format!("{:.6}", density[cp.position as usize]);
                let laplacian_text = format!(
                    "{:.6}",
                    laplacian(cp.position as usize, density, grid)
                );
                // each atom of the point gets a row, once per atom even if
                // the point involves more than one image of it
                for (i, atom) in cp.atoms.iter().enumerate() {
                    if cp.atoms[..i]
                        .iter()
                        .any(|a| a.atom_index() == atom.atom_index())
                    {
                        continue;
                    }
                    let members = cp
                        .atoms
                        .iter()
                        .map(|a| a.image_sub(atom.image()))
                        .filter(|a| {
                            *a != EncodedAtom::new_zero_image(atom.atom_index())
                        })
                        .map(format_member)
                        .collect::<Vec<String>>()
                        .join(", ");
                    let atom = atom.atom_index() as usize;
                    if kind == "Bond" {
                        coordination[atom] += 1;
                    }
                    rows[atom].push(CriticalPointRow {
                        kind,
                        position: position.clone(),
                        density: density_text.clone(),
                        laplacian: laplacian_text.clone(),
                        members,
                    });
                }
            }
        }
        let positions = atoms
            .positions
            .iter()
            .map(|p| file_type.coordinate_format(*p))
            .collect();
        let summary = format!(
            "## Topological Information\n\n * Total nuclei: {}\n * Total bonds: {}\n * Total rings: {}\n * Total cages: {}\n",
            critical_points.nuclei.len(),
            critical_points.bonds.len(),
            critical_points.rings.len(),
            critical_points.cages.len()
        );
        Self {
            summary,
            positions,
            coordination,
            rows,
        }
    }
}

/// Formats an atom as its 1-based index followed by its image if not `[0, 0, 0]`.
fn format_member(atom: EncodedAtom) -> String {
    let index = atom.atom_index() + 1;
    if atom.image().is_zero() {
        format!("{}", index)
    } else {
        let [x, y, z] = atom.image().decode();
        format!("{}({} {} {})", index, x, y, z)
    }
}

impl fmt::Display for CriticalPointOutput {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.summary)?;
        for (i, (position, rows)) in
            self.positions.iter().zip(&self.rows).enumerate()
        {
            writeln!(
                f,
                "\n### Atom: {}\n * Position: {:.6} {:.6} {:.6}\n * Coordination: {}\n * Critical Points:",
                i + 1,
                position[0],
                position[1],
                position[2],
                self.coordination[i],
            )?;
            // widths are the longest entry in each column or its title
            let width =
                |title: &str, column: &dyn Fn(&CriticalPointRow) -> usize| {
                    rows.iter().map(column).max().unwrap_or(0).max(title.len())
                };
            let kind = width("Type", &|r| r.kind.len());
            let pos = width("", &|r| {
                r.position.iter().map(String::len).max().unwrap()
            });
            let den = width("Density", &|r| r.density.len());
            let lap = width("Laplacian", &|r| r.laplacian.len());
            let mem = width("Members", &|r| r.members.len());
            let pos_title = (pos * 3 + 2).max("Position".len());
            writeln!(
                f,
                "| {:^kind$} | {:^pos_title$} | {:^den$} | {:^lap$} | {:^mem$} |",
                "Type", "Position", "Density", "Laplacian", "Members"
            )?;
            writeln!(
                f,
                "|-{:-^kind$}-|-{:-^pos_title$}-|-{:-^den$}-|-{:-^lap$}-|-{:-^mem$}-|",
                "", "", "", "", ""
            )?;
            for r in rows {
                // right align the coordinates within the position column
                let position = format!(
                    "{:>pos$} {:>pos$} {:>pos$}",
                    r.position[0], r.position[1], r.position[2]
                );
                writeln!(
                    f,
                    "| {:^kind$} | {:>pos_title$} | {:>den$} | {:>lap$} | {:<mem$} |",
                    r.kind, position, r.density, r.laplacian, r.members
                )?;
            }
        }
        Ok(())
    }
}

/// The partitioned charge of each atom formatted as a markdown table by its
/// [`Display`](fmt::Display) implementation.
pub struct PartitionTable {
    /// The column titles.
    titles: Vec<String>,
    /// The rows of the table as a vector of strings.
    rows: Vec<Vec<String>>,
    /// The summary of the vacuum and totals.
    footer: String,
}

impl PartitionTable {
    /// Creates the table from the partition.
    ///
    /// # Arguments
    /// * `partition`: The integrated Bader volumes.
    /// * `boundary_voxels`: The number of voxels shared between volumes.
    /// * `total_voxels`: The number of voxels in the grid.
    pub fn new(
        partition: &Partition,
        boundary_voxels: usize,
        total_voxels: usize,
    ) -> Self {
        let atoms = partition.atoms();
        let names = density_names(partition.charge[0].len());
        let titles = ["Atom #"]
            .iter()
            .chain(names)
            .chain(&["Volume", "Distance", "Error"])
            .map(|s| s.to_string())
            .collect();
        let rows = (0..atoms)
            .map(|i| {
                [format!("{}", i + 1)]
                    .into_iter()
                    .chain(
                        partition.charge[i].iter().map(|c| format!("{:.6}", c)),
                    )
                    .chain([
                        format!("{:.6}", partition.volume[i]),
                        format!("{:.6}", partition.radius[i]),
                        format!("{:.6}", partition.error[i]),
                    ])
                    .collect()
            })
            .collect();
        // the last entries are the vacuum
        let vacuum_charge = &partition.charge[atoms];
        let vacuum_volume = partition.volume[atoms];
        let total_charge = partition.charge[..atoms].iter().fold(
            vec![0.0; names.len()],
            |mut total, charge| {
                total.iter_mut().zip(charge).for_each(|(t, c)| *t += c);
                total
            },
        );
        let total_volume: f64 = partition.volume[..atoms].iter().sum();
        let mut footer = String::new();
        let mut line = |label: String, value: String| {
            footer.push_str(&format!("\n * {:<20}{:>13}", label + ":", value));
        };
        for (name, charge) in names.iter().zip(vacuum_charge) {
            line(format!("Vacuum {}", name), format!("{:.4}", charge));
        }
        line(
            String::from("Vacuum Volume"),
            format!("{:.4}", vacuum_volume),
        );
        for (name, charge) in names.iter().zip(&total_charge) {
            line(format!("Partitioned {}", name), format!("{:.4}", charge));
        }
        line(
            String::from("Partitioned Volume"),
            format!("{:.4}", total_volume),
        );
        line(
            String::from("Boundary Voxels"),
            format!("{}", boundary_voxels),
        );
        line(String::from("Total Voxels"), format!("{}", total_voxels));
        Self {
            titles,
            rows,
            footer,
        }
    }
}

impl fmt::Display for PartitionTable {
    /// Creates a String representation of the Table.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let widths = self
            .titles
            .iter()
            .enumerate()
            .map(|(i, title)| {
                self.rows
                    .iter()
                    .map(|row| row[i].len())
                    .max()
                    .unwrap_or(0)
                    .max(title.len())
            })
            .collect::<Vec<usize>>();
        writeln!(f, "## Partition Information")?;
        write!(f, "|")?;
        for (title, width) in self.titles.iter().zip(&widths) {
            write!(f, " {:^width$} |", title)?;
        }
        write!(f, "\n|")?;
        for width in &widths {
            write!(f, "-{:-^width$}-|", "")?;
        }
        for row in &self.rows {
            write!(f, "\n|")?;
            for (cell, width) in row.iter().zip(&widths) {
                write!(f, " {:>width$} |", cell)?;
            }
        }
        write!(f, "\n{}", self.footer)
    }
}

/// Write a string to filename. Creates a new file regardless of what exists.
pub fn write(string: String, filename: String) -> std::io::Result<()> {
    File::create(filename)?.write_all(string.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::Lattice;
    use crate::critical::CriticalPointKind;
    use crate::io::cube::Cube;
    use crate::voxel_map::EncodedImage;

    #[test]
    fn test_density_names() {
        assert_eq!(density_names(1), &["Charge"]);
        assert_eq!(density_names(2), &["Charge", "Spin"]);
        assert_eq!(density_names(4), &["Charge", "Spin X", "Spin Y", "Spin Z"]);
    }

    #[test]
    fn test_partition_table() {
        let partition = Partition {
            charge: vec![vec![1.5, 0.25], vec![2.5, -0.25], vec![0.1, 0.0]],
            volume: vec![10.0, 20.0, 0.5],
            radius: vec![1.0, 2.0],
            error: vec![0.001, -0.001, 0.0],
        };
        let table = PartitionTable::new(&partition, 12, 100).to_string();
        let expected = "\
## Partition Information
| Atom # |  Charge  |   Spin    |  Volume   | Distance |   Error   |
|--------|----------|-----------|-----------|----------|-----------|
|      1 | 1.500000 |  0.250000 | 10.000000 | 1.000000 |  0.001000 |
|      2 | 2.500000 | -0.250000 | 20.000000 | 2.000000 | -0.001000 |

 * Vacuum Charge:             0.1000
 * Vacuum Spin:               0.0000
 * Vacuum Volume:             0.5000
 * Partitioned Charge:        4.0000
 * Partitioned Spin:          0.0000
 * Partitioned Volume:       30.0000
 * Boundary Voxels:               12
 * Total Voxels:                 100";
        assert_eq!(table, expected);
    }

    #[test]
    fn test_critical_point_output() {
        let lattice = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        let atoms = Atoms::new(
            Lattice::new(lattice),
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            String::new(),
        );
        let grid = Grid::new([4, 4, 4], lattice, [0.0; 3]);
        let density = vec![1.0; 64];
        let nuclei = vec![
            CriticalPoint::new(
                0,
                CriticalPointKind::Nuclei,
                vec![EncodedAtom::new_zero_image(0)].into(),
            ),
            CriticalPoint::new(
                32,
                CriticalPointKind::Nuclei,
                vec![EncodedAtom::new_zero_image(1)].into(),
            ),
        ];
        // one bond inside the cell and one across the periodic boundary
        let bonds = vec![
            CriticalPoint::new(
                16,
                CriticalPointKind::Bond,
                vec![
                    EncodedAtom::new_zero_image(0),
                    EncodedAtom::new_zero_image(1),
                ]
                .into(),
            ),
            CriticalPoint::new(
                48,
                CriticalPointKind::Bond,
                vec![
                    EncodedAtom::new_zero_image(1),
                    EncodedAtom::new(0, EncodedImage::new([1, 0, 0])),
                ]
                .into(),
            ),
        ];
        let critical_points = CriticalPoints {
            nuclei,
            bonds,
            ..Default::default()
        };
        let output = CriticalPointOutput::new(
            &atoms,
            &critical_points,
            &density,
            &grid,
            FileType::Cube(Cube {}),
        )
        .to_string();
        let expected = "\
## Topological Information

 * Total nuclei: 2
 * Total bonds: 2
 * Total rings: 0
 * Total cages: 0

### Atom: 1
 * Position: 0.000000 0.000000 0.000000
 * Coordination: 2
 * Critical Points:
|  Type   |          Position          | Density  | Laplacian |  Members  |
|---------|----------------------------|----------|-----------|-----------|
| Nucleus | 0.000000 0.000000 0.000000 | 1.000000 |  0.000000 |           |
|  Bond   | 1.000000 0.000000 0.000000 | 1.000000 |  0.000000 | 2         |
|  Bond   | 3.000000 0.000000 0.000000 | 1.000000 |  0.000000 | 2(-1 0 0) |

### Atom: 2
 * Position: 2.000000 0.000000 0.000000
 * Coordination: 2
 * Critical Points:
|  Type   |          Position          | Density  | Laplacian | Members  |
|---------|----------------------------|----------|-----------|----------|
| Nucleus | 2.000000 0.000000 0.000000 | 1.000000 |  0.000000 |          |
|  Bond   | 1.000000 0.000000 0.000000 | 1.000000 |  0.000000 | 1        |
|  Bond   | 3.000000 0.000000 0.000000 | 1.000000 |  0.000000 | 1(1 0 0) |
";
        assert_eq!(output, expected);
    }

    #[test]
    fn test_format_member() {
        assert_eq!(format_member(EncodedAtom::new_zero_image(4)), "5");
        assert_eq!(
            format_member(EncodedAtom::new(4, EncodedImage::new([0, -1, 2]))),
            "5(0 -1 2)"
        );
    }
}
