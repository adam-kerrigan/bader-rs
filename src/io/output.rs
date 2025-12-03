use crate::critical::CriticalPoint;
use std::fs::File;
use std::io::Write;

pub struct OutputFile {
    partition_table: PartitionTable,
    critical_points: CriticalPointOutput,
}

impl OutputFile {
    pub fn new(
        partition_table: PartitionTable,
        critical_points: CriticalPointOutput,
    ) -> Self {
        Self {
            partition_table,
            critical_points,
        }
    }
}

struct CriticalPointOutput {
    ordered_postitions: Vec<String>,
    atom_nuclei_map: Vec<Vec<CriticalPoint>>,
    atom_bond_map: Vec<Vec<CriticalPoint>>,
    atom_ring_map: Vec<Vec<CriticalPoint>>,
    atom_cage_map: Vec<Vec<CriticalPoint>>,
}

/// Structure that contains and builds the table.
pub struct PartitionTable {
    /// How wide each column is.
    column_width: Vec<usize>,
    /// The number of charge and spin densities.
    density_num: usize,
    /// The rows of the table as a vector of strings.
    rows: Vec<Vec<String>>,
    footer: String,
}

impl PartitionTable {
    /// Creates a new structure and sets the minimum widths of each.
    pub fn new(
        partitioned_density: &[Box<[f64]>],
        partitioned_volume: &[f64],
        radius: &[f64],
        errors: &[f64],
        boundary_voxels: usize,
        total_voxels: usize,
    ) -> Self {
        let density_num = partitioned_density[0].len();
        let mut rows: Vec<Vec<String>> =
            vec![Vec::with_capacity(4 + density_num)];
        let mut column_width: Vec<usize> = vec![6; 4 + density_num];
        column_width[2 + density_num] = 8;
        // calculate the total density for each density supplied
        let total_density: Vec<f64> = partitioned_density.iter().fold(
            vec![0.0; partitioned_density[0].len()],
            |mut sum, d| {
                sum.iter_mut().zip(d).for_each(|(tpd, pd)| *tpd += pd);
                sum
            },
        );
        // the last value is is the vacuum and it has definitely been added
        let vacuum_density = partitioned_density.last().unwrap();
        let total_partitioned_density = total_density
            .iter()
            .zip(vacuum_density)
            .map(|(td, vd)| td - vd)
            .collect::<Vec<f64>>();
        // the volume is the same for all densities
        let total_volume: f64 = partitioned_volume.iter().sum();
        // the last value is is the vacuum and it has definitely been added
        let vacuum_volume = *partitioned_volume.last().unwrap();
        let total_partitioned_volume = total_volume - vacuum_volume;
        partitioned_density
            .iter()
            .zip(partitioned_volume)
            .zip(radius)
            .zip(errors)
            .for_each(|(((density, volume), radius), error)| {
                let mut row: Vec<String> =
                    Vec::with_capacity(rows[0].capacity());
                row.push(format!("{}", rows.len() - 1));
                density.iter().for_each(|d| row.push(format!("{:.6}", d)));
                row.push(format!("{:.6}", volume));
                row.push(format!("{:.6}", radius));
                row.push(format!("{:.6}", error));
                for (i, col) in row.iter().enumerate() {
                    column_width[i] = column_width[i].max(col.len());
                }
                rows.push(row);
            });
        let mut footer = match density_num.cmp(&2) {
            std::cmp::Ordering::Less => format!(
                "\n * Vacuum Charge: {:>18.4}\n * Vacuum Volume: {:>18.4}\n * Partitioned Charge: {:>13.4}\n * Partitioned Volume: {:>13.4}",
                vacuum_density[0],
                vacuum_volume,
                total_partitioned_density[0],
                total_partitioned_volume,
            ),
            std::cmp::Ordering::Equal => format!(
                "\n * Vacuum Charge: {:>18.4}\n * Vacuum Spin: {:>20.4}\n * Vacuum Volume: {:>18.4}\n * Partitioned Charge: {:>13.4}\n * Partitioned Spin: {:>15.4}\n * Partitioned Volume: {:>13.4}",
                vacuum_density[0],
                vacuum_density[1],
                vacuum_volume,
                total_partitioned_density[0],
                total_partitioned_density[1],
                total_partitioned_volume,
            ),
            std::cmp::Ordering::Greater => format!(
                "\n * Vacuum Charge: {:>18.4}\n * Vacuum Spin X: {:>18.4}\n * Vacuum Spin Y: {:>18.4}\n * Vacuum Spin Z: {:>18.4}\n * Vacuum Volume: {:>18.4}\n * Partitioned Charge: {:>13.4}\n * Partitioned Spin X: {:>13.4}\n * Partitioned Spin Y: {:>13.4}\n * Partitioned Spin Z: {:>13.4}\n * Partitioned Volume: {:>13.4}",
                vacuum_density[0],
                vacuum_density[1],
                vacuum_density[2],
                vacuum_density[3],
                vacuum_volume,
                total_partitioned_density[0],
                total_partitioned_density[1],
                total_partitioned_density[2],
                total_partitioned_density[3],
                total_partitioned_volume,
            ),
        };
        footer.push_str(&format!(
            "\n * Boundary Voxels: {:>16.4}\n * Total Voxels: {:>19.4}",
            boundary_voxels, total_voxels
        ));
        Self {
            column_width,
            density_num,
            rows,
            footer,
        }
    }

    /// Creates and formats the header.
    fn format_header(&self) -> String {
        let mut header = String::new();
        let mut iter = self.column_width.iter();
        header.push_str(&format!(
            "| {:^width$} |",
            "Atom #",
            width = iter.next().unwrap()
        ));
        header.push_str(&format!(
            " {:^width$} |",
            "Charge",
            width = iter.next().unwrap()
        ));
        match self.density_num.cmp(&2) {
            std::cmp::Ordering::Equal => {
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin",
                    width = iter.next().unwrap()
                ));
            }
            std::cmp::Ordering::Greater => {
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin X",
                    width = iter.next().unwrap()
                ));
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin Y",
                    width = iter.next().unwrap()
                ));
                header.push_str(&format!(
                    " {:^width$} |",
                    "Spin Z",
                    width = iter.next().unwrap()
                ));
            }
            std::cmp::Ordering::Less => (),
        }
        header.push_str(&format!(
            " {:^width$} |",
            "Volume",
            width = iter.next().unwrap()
        ));
        header.push_str(&format!(
            " {:^width$} |",
            "Distance",
            width = iter.next().unwrap()
        ));
        header.push_str(&format!(
            " {:^width$} |\n",
            "Error",
            width = iter.next().unwrap()
        ));
        header
    }

    /// Creates and formats a separator.
    fn format_separator(&self) -> String {
        let mut separator = String::from("|");
        self.column_width.iter().for_each(|w| {
            separator.push_str(&format!("-{:-^width$}-|", "-", width = w));
        });
        separator
    }

    /// Creates a String representation of the Table.
    pub fn get_string(self) -> String {
        let mut table = String::new();
        table.push_str(&self.format_header());
        self.rows.iter().for_each(|r| {
            if r.is_empty() {
                table.push_str(&self.format_separator());
            } else {
                let mut row = String::from("|");
                r.iter().zip(&self.column_width).for_each(|(s, w)| {
                    row.push_str(&format!(" {:>width$} |", s, width = w))
                });
                table.push_str(&row);
            }
            table.push('\n');
        });
        table.push_str(&self.footer);
        table
    }
}

/// Write a string to filename. Creates a new file regardless of what exists.
pub fn write(string: String, filename: String) -> std::io::Result<()> {
    let mut bader_file = File::create(filename)?;
    bader_file.write_all(string.as_bytes())?;
    Ok(())
}
