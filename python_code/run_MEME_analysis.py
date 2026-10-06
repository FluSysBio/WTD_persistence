"""Script to run HyPhy MEME analyses for Alpha and Delta gene trees."""

import argparse
import csv
import logging
import re
import subprocess
import sys
from pathlib import Path

P_VALUE = "0.1"

SUMMARY_COLUMNS = [
    "gene", "clusters", "singletons", "deer_branches",
    "bg_branches", "fg_length", "bg_length", "bg_pct",
]


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def get_analysis_name(directory_name):
    """Extract variant and gene from a directory name."""
    name = directory_name.removesuffix("-gene")
    parts = name.split("_", 1)

    if len(parts) != 2:
        return None

    variant, gene = parts
    if variant.lower() not in {"alpha", "delta"}:
        return None

    return f"{variant.capitalize()}_{gene}"


def find_input_files(gene_dir):
    """Find the alignment and phylogenetic tree."""
    alignments = sorted(gene_dir.glob("*.fasta"))
    trees = sorted(gene_dir.glob("*.treefile"))

    if not alignments or not trees:
        return None

    return alignments[0], trees[0]


def run_command(command, log_file):
    """Execute a command and save its output to a log file."""
    with log_file.open("w") as log:
        result = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )

    if result.returncode != 0:
        logging.error("Command failed. See log: %s", log_file)
        return False

    return True


def prepare_clusters(tree, threshold, scripts_dir, output_dir, name):
    """Label foreground branches and identify clusters."""
    labeled_tree = output_dir / f"{name}_labeled.nwk"
    clusters_file = output_dir / f"{name}_clusters.tsv"

    command = [
        sys.executable,
        str(scripts_dir / "prepare_clusters.py"),
        "--tree", str(tree),
        "--threshold", str(threshold),
        "--out-tree", str(labeled_tree),
        "--out-clusters", str(clusters_file),
    ]

    success = run_command(
        command, output_dir / f"{name}_clusters.log"
    )

    if not success or not labeled_tree.is_file() or not clusters_file.is_file():
        logging.error("Cluster preparation failed for %s", name)
        return None

    return labeled_tree, clusters_file


def summarize_branches(labeled_tree, clusters_file, name):
    """Summarize foreground/background branch counts and lengths."""
    tree_text = labeled_tree.read_text()
    branches = re.findall(r"(\{Deer\})?:([0-9.eE+-]+)", tree_text)

    foreground = [float(length) for tag, length in branches if tag]
    background = [float(length) for tag, length in branches if not tag]

    fg_length = sum(foreground)
    bg_length = sum(background)
    total_length = fg_length + bg_length

    with clusters_file.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))[1:]

    cluster_ids = [row[0] for row in rows if row]

    return {
        "gene": name,
        "clusters": len({
            cluster for cluster in cluster_ids if cluster != "singleton"
        }),
        "singletons": sum(
            cluster == "singleton" for cluster in cluster_ids
        ),
        "deer_branches": len(foreground),
        "bg_branches": len(background),
        "fg_length": f"{fg_length:.5f}",
        "bg_length": f"{bg_length:.5f}",
        "bg_pct": f"{100 * bg_length / total_length:.1f}"
        if total_length else "0.0",
    }


def run_meme(alignment, tree, output_file, log_file, foreground=False):
    """Run MEME on all branches or Deer foreground branches."""
    command = [
        "hyphy", "meme",
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--pvalue", P_VALUE,
        "--output", str(output_file),
    ]

    if foreground:
        command.extend(["--branches", "Deer"])

    success = run_command(command, log_file)

    if success and not output_file.is_file():
        logging.error("MEME output not found: %s", output_file)
        return False

    return success


def process_gene(gene_dir, output_dir, scripts_dir, threshold):
    """Run cluster preparation and MEME analyses for one gene."""
    name = get_analysis_name(gene_dir.name)

    if not name:
        logging.warning("Skipping unrecognized directory: %s", gene_dir.name)
        return None

    inputs = find_input_files(gene_dir)
    if inputs is None:
        logging.warning("Missing FASTA or tree file: %s", gene_dir)
        return None

    alignment, tree = inputs
    logging.info("Processing %s", name)

    cluster_outputs = prepare_clusters(
        tree, threshold, scripts_dir, output_dir, name
    )
    if cluster_outputs is None:
        return None

    labeled_tree, clusters_file = cluster_outputs
    summary = summarize_branches(labeled_tree, clusters_file, name)

    run_meme(
        alignment, tree,
        output_dir / f"{name}_MEME_all.json",
        output_dir / f"{name}_MEME_all.log",
    )

    if "{Deer}" in labeled_tree.read_text():
        run_meme(
            alignment, labeled_tree,
            output_dir / f"{name}_MEME_deer.json",
            output_dir / f"{name}_MEME_deer.log",
            foreground=True,
        )
    else:
        logging.info("No Deer branches for %s; skipping foreground MEME", name)

    logging.info("Finished %s", name)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Run HyPhy MEME analyses for Alpha and Delta genes."
    )
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("--threshold", type=float, default=0.0015)
    parser.add_argument(
        "--scripts",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing prepare_clusters.py",
    )
    args = parser.parse_args()

    setup_logging()

    data_dir = args.data_dir.resolve()
    scripts_dir = args.scripts.resolve()

    if not data_dir.is_dir():
        parser.error(f"Directory not found: {data_dir}")

    output_dir = data_dir / "meme_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_dirs = sorted(
        path for path in data_dir.glob("*-gene") if path.is_dir()
    )
    if not gene_dirs:
        parser.error(f"No *-gene directories found in {data_dir}")

    summary_rows = []

    for gene_dir in gene_dirs:
        row = process_gene(
            gene_dir, output_dir, scripts_dir, args.threshold
        )
        if row is not None:
            summary_rows.append(row)

    summary_file = output_dir / "branch_sets_summary.tsv"

    with summary_file.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    logging.info(
        "Processed %d of %d directories", len(summary_rows), len(gene_dirs)
    )
    logging.info("Results saved to %s", output_dir)
    print(summary_file.read_text())


if __name__ == "__main__":
    main()