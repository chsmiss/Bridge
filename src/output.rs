use crate::assembler::AssemblyProduct;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write_outputs(product: &AssemblyProduct, output_dir: &Path) -> Result<()> {
    write_fasta(
        &output_dir.join("primary_contigs.fasta"),
        product
            .primary_sequences
            .iter()
            .enumerate()
            .map(|(index, sequence)| {
                (
                    format!("contig_{:06} len={}", index + 1, sequence.len()),
                    sequence.as_slice(),
                )
            }),
    )?;

    write_fasta(
        &output_dir.join("unitigs.fasta"),
        product.unitig_graph.unitigs.iter().map(|unitig| {
            (
                format!(
                    "unitig_{:06} len={} cov={:.3} min_cov={} max_cov={} circular={}",
                    unitig.id + 1,
                    unitig.length,
                    unitig.mean_coverage,
                    unitig.min_coverage,
                    unitig.max_coverage,
                    unitig.circular
                ),
                unitig.sequence.as_slice(),
            )
        }),
    )?;

    write_gfa(product, &output_dir.join("assembly.gfa"))?;
    let stats_file = File::create(output_dir.join("run_profile.json"))
        .context("failed to create run_profile.json")?;
    serde_json::to_writer_pretty(BufWriter::new(stats_file), &product.stats)
        .context("failed to write run_profile.json")?;
    Ok(())
}

fn write_fasta<'a, I>(path: &Path, records: I) -> Result<()>
where
    I: IntoIterator<Item = (String, &'a [u8])>,
{
    let file = File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for (header, sequence) in records {
        writeln!(writer, ">{header}")?;
        for chunk in sequence.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_gfa(product: &AssemblyProduct, path: &Path) -> Result<()> {
    let file = File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "H\tVN:Z:1.0")?;
    for unitig in &product.unitig_graph.unitigs {
        writeln!(
            writer,
            "S\tu{}\t{}\tLN:i:{}\tKC:f:{:.3}",
            unitig.id,
            String::from_utf8_lossy(&unitig.sequence),
            unitig.length,
            unitig.mean_coverage
        )?;
    }
    let mut transitions: Vec<_> = product.transitions.iter().collect();
    transitions.sort_unstable_by_key(|(edge, _)| **edge);
    for (&(source, target), evidence) in transitions {
        if evidence.direct_reads == 0 {
            continue;
        }
        writeln!(
            writer,
            "L\tu{}\t+\tu{}\t+\t{}M\tDR:i:{}\tPE:i:{}",
            source,
            target,
            product.unitig_graph.k,
            evidence.direct_reads,
            evidence.read_pairs
        )?;
    }
    writer.flush()?;
    Ok(())
}
