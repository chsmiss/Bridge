from pathlib import Path
import re

assembler = Path("src/assembler.rs")
text = assembler.read_text()
pattern = re.compile(
    r"    let mut excluded = FxHashSet::default\(\);\n"
    r"    for alleles in groups\.values\(\) \{\n"
    r"        let Some\(primary\) = alleles\.iter\(\)\.copied\(\)\.max_by\(\|left, right\| \{.*?"
    r"        for allele in alleles \{\n"
    r"            if allele\.unitig_id != primary\.unitig_id \{\n"
    r"                excluded\.insert\(allele\.unitig_id\);\n"
    r"            \}\n"
    r"        \}\n"
    r"    \}\n",
    re.S,
)
replacement = '''    let mut excluded = FxHashSet::default();
    for alleles in groups.values() {
        // Sharing graph boundaries is not enough to prove a biological
        // bubble. Collapse alternatives only when at least two alleles have
        // independent direct-read support on both flanks. This prevents
        // repeat/orientation artefacts from truncating a linear primary path.
        let supported: Vec<&BubbleAllele> = alleles
            .iter()
            .copied()
            .filter(|allele| allele.physically_flanked)
            .collect();
        if supported.len() < 2 {
            continue;
        }
        let Some(primary) = supported.iter().copied().max_by(|left, right| {
            left.mean_coverage
                .total_cmp(&right.mean_coverage)
                .then_with(|| {
                    (left.left_support + left.right_support)
                        .cmp(&(right.left_support + right.right_support))
                })
                .then_with(|| right.unitig_id.cmp(&left.unitig_id))
        }) else {
            continue;
        };
        for allele in supported {
            if allele.unitig_id != primary.unitig_id {
                excluded.insert(allele.unitig_id);
            }
        }
    }
'''
text, count = pattern.subn(replacement, text, count=1)
if count != 1:
    raise SystemExit(f"bubble block replacements: {count}")
assembler.write_text(text)

# These files are temporary implementation scaffolding and must not remain in
# the production branch. The workflow itself is removed separately through
# the GitHub API because Actions tokens cannot modify workflow files.
Path(".github/patches/bubble-safety.patch").unlink(missing_ok=True)
Path(".github/scripts/apply_bubble_fix.py").unlink(missing_ok=True)
