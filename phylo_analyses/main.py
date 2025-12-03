#!/usr/bin/env python3
"""
  python main.py -i "D:\source\thesis\test_mafft\output2.fasta"
"""

import argparse
import os
import subprocess
import shutil
import sys
from Bio import Phylo

# ---------- CONFIG ----------
IQTREE_BIN_DEFAULT = "iqtree2"
UFB = 1000
ALRT = 1000
THREADS = "AUTO"
# -----------------------------

def which_or_err(bin_name, hint=None):
    path = shutil.which(bin_name)
    if path:
        return path
    raise FileNotFoundError(f"'{bin_name}' is not found in PATH. {hint or ''}")

def run_iqtree(alignment_file, iqtree_bin, ufb, alrt, threads, prefix):
    iqtree_path = which_or_err(iqtree_bin, "No IQ-TREE2.")
    cmd = [
        iqtree_path,
        "-s", alignment_file,
        "-m", "MFP",
        "-bb", str(ufb),
        "-alrt", str(alrt),
        "-nt", str(threads),
        "-pre", prefix, 
        "-redo"      
    ]
    print("[IQ-TREE] Uruchamiam:", " ".join(cmd))
    workdir = os.path.dirname(os.path.abspath(alignment_file))
    proc = subprocess.run(cmd, cwd=workdir)
    print("RC:", proc.returncode)
    if proc.returncode != 0:
        raise RuntimeError("IQ-TREE error.")


def parse_best_model(iqtree_report):
    if not os.path.exists(iqtree_report):
        raise FileNotFoundError(f"No file .iqtree: {iqtree_report}")

    best = None
    with open(iqtree_report) as fh:
        for line in fh:
            if "Best-fit model" in line or "Best model" in line:
                parts = line.split(":")
                if len(parts) >= 2:
                    model = parts[1].strip().split()[0]
                    best = model
                    break

    if not best:
        raise RuntimeError("Best model not found in .iqtree file.")
    print("[MODEL] Best model:", best)
    return best

def find_treefile(prefix):
    for ext in [".contree", ".treefile", ".treefile.gz"]:
        f = prefix + ext
        if os.path.exists(f):
            return f
    raise FileNotFoundError("Files not found .contree/.treefile")

def plot_static_tree(treefile, out_png):
    print("[PLOT] static tree => PNG:", out_png)
    tree = Phylo.read(treefile, "newick")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(1,1,1)
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print("[PLOT] PNG saved.")

def open_ete3(treefile):
    try:
        from ete3 import Tree, TreeStyle, NodeStyle
    except Exception:
        print("[ETE3] No ETE3 => skipping viewer.")
        return

<<<<<<< HEAD
# Save the tree in Newick format
Phylo.write(tree, "Phylogenetic_Tree_Construction.nwk", "nexml")  
=======
    print("[ETE3] Opening viewer.")
    try:
        t = Tree(treefile, format=1)
    except:
        t = Tree(treefile)
>>>>>>> model

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = True

    for n in t.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        n.set_style(nstyle)

    t.show(tree_style=ts)

def main():
    p = argparse.ArgumentParser(description="IQ-TREE pipeline - aligment prepared.")
    p.add_argument("-i", "--input", required=True, help="Prepared Alignment FASTA.")
    p.add_argument("--iqtree-bin", default=IQTREE_BIN_DEFAULT)
    p.add_argument("--threads", default=THREADS)
    p.add_argument("--no-viewer", action="store_true")
    args = p.parse_args()

    alignment = args.input
    if not os.path.exists(alignment):
        print("No Aligment file:", alignment)
        sys.exit(1)

    alignment = os.path.abspath(alignment)

    prefix = os.path.splitext(alignment)[0]

    run_iqtree(alignment, args.iqtree_bin, UFB, ALRT, args.threads, prefix)

    model = parse_best_model(prefix + ".iqtree")
    treefile = find_treefile(prefix)
    print("[TREE] tree file:", treefile)

    png = prefix + ".png"
    plot_static_tree(treefile, png)


    if not args.no_viewer:
        open_ete3(treefile)

if __name__ == "__main__":
    main()



##################################################### warpper mafft ############################################################################



# def run_mafft(input_fasta, output_alignment, mafft_bin=MAFFT_DEFAULT):
#     mafft_path = which_or_err(mafft_bin, "nstall mafft")
#     cmd = [mafft_path, "--auto", input_fasta]
#     print("[MAFFT] Uruchamiam MAFFT:", " ".join(cmd))
#     with open(output_alignment, "w") as out:
#         proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
#     if proc.returncode != 0:
#         print("[MAFFT] stderr:", proc.stderr)
#         raise RuntimeError("MAFFT error.")
#     print(f"[MAFFT] Aligment fil saved: {output_alignment}")
