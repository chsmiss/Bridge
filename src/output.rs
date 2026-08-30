use crate::assembler::{AssemblyProduct, BubbleAllele};
use crate::dna::reverse_complement;
use crate::scaffold::ScaffoldLink;
use anyhow::{Context, Result};
use rustc_hash::FxHashSet;
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
        &output_dir.join("primary_scaffolds.fasta"),
        product
            .scaffold_sequences
            .iter()
            .enumerate()
            .map(|(index, sequence)| {
                (
                    format!("scaffold_{:06} len={}", index + 1, sequence.len()),
                    sequence.as_slice(),
                )
            }),
    )?;
    write_scaffold_links(
        &output_dir.join("scaffold_links.tsv"),
        &product.scaffold_links,
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

    write_bubble_fasta(
        &output_dir.join("variants.fasta"),
        &product.bubble_alleles,
        false,
    )?;
    write_bubble_fasta(
        &output_dir.join("haplotigs.fasta"),
        &product.bubble_alleles,
        true,
    )?;
    write_bubble_table(
        &output_dir.join("bubble_alleles.tsv"),
        &product.bubble_alleles,
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
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
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

fn write_bubble_fasta(path: &Path, alleles: &[BubbleAllele], haplotigs: bool) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    let mut output_index = 0_usize;

    for allele in alleles {
        let sequence = if haplotigs {
            let Some(sequence) = allele.haplotig_sequence.as_deref() else {
                continue;
            };
            sequence
        } else {
            allele.allele_sequence.as_slice()
        };
        let reverse = reverse_complement(sequence);
        let canonical = if reverse.as_slice() < sequence {
            reverse
        } else {
            sequence.to_vec()
        };
        if !seen.insert(canonical.clone()) {
            continue;
        }
        output_index += 1;
        writeln!(
            writer,
            ">{}_{:06} bubble={} allele={} len={} cov={:.3} left_reads={} right_reads={} flanked={}",
            if haplotigs { "haplotig" } else { "variant" },
            output_index,
            allele.bubble_id,
            allele.allele_index,
            canonical.len(),
            allele.mean_coverage,
            allele.left_support,
            allele.right_support,
            allele.physically_flanked
        )?;
        for chunk in canonical.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_scaffold_links(path: &Path, links: &[ScaffoldLink]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "source_component\ttarget_component\tpair_support\tgap_bases"
    )?;
    for link in links {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            link.source_component, link.target_component, link.pair_support, link.gap_bases
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_bubble_table(path: &Path, alleles: &[BubbleAllele]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "bubble_id\tallele_index\tunitig_id\tlength\tmean_coverage\tleft_support\tright_support\tphysically_flanked\tpath"
    )?;
    for allele in alleles {
        let path_text = allele
            .path
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}",
            allele.bubble_id,
            allele.allele_index,
            allele.unitig_id,
            allele.length,
            allele.mean_coverage,
            allele.left_support,
            allele.right_support,
            allele.physically_flanked,
            path_text
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_gfa(product: &AssemblyProduct, path: &Path) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
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
    for source in 0..product.unitig_graph.unitigs.len() as u32 {
        for edge_index in product.unitig_graph.out_range(source) {
            let target = product.unitig_graph.out_targets[edge_index];
            let evidence = product
                .transitions
                .get(&(source, target))
                .copied()
                .unwrap_or_default();
            writeln!(
                writer,
                "L\tu{}\t+\tu{}\t+\t{}M\tDR:i:{}\tGR:i:{}\tPE:i:{}",
                source,
                target,
                product.unitig_graph.k,
                evidence.direct_reads,
                evidence.gapped_reads,
                evidence.read_pairs
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}
