
from email import header
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
import customtkinter as ctk
import os
import subprocess
from Bio import Phylo
import matplotlib.pyplot as plt
from pathlib import Path
import threading
from datetime import datetime
import platform

ctk.set_appearance_mode("light")
ctk.set_default_color_theme("blue")

class PhylogenyApp:
    def __init__(self):
        self.root = ctk.CTk()
        self.root.title("Phylogenetic Analysis Tool")
        self.root.geometry("800x600")
        
        self.fasta_files = []
        self.base_output_dir = Path.home() / "phylogeny_analysis"
        self.base_output_dir.mkdir(exist_ok=True)

        self.input_for_mafft = None
        self.study_dir = None
        self.dir_input = None
        self.dir_alignment = None
        self.dir_tree = None
        self.dir_iqtree = None
        self.dir_plots = None

        self.setup_ui()
    
    def create_study_dirs(self):
        assert self.base_output_dir.exists()

        dlg = ctk.CTkInputDialog(
            text="Enter a name for this study:",
            title="Study Name"
        )
        study_name = dlg.get_input()
        if not study_name:
            raise RuntimeError("No study name provided")

        safe_name = study_name.strip().replace(" ", "_")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.study_dir = self.base_output_dir / f"{safe_name}_{timestamp}"
        
        self.dir_input = self.study_dir / "input"
        self.dir_alignment = self.study_dir / "alignment"
        self.dir_tree = self.study_dir / "tree"
        self.dir_iqtree = self.dir_tree / "iqtree"
        self.dir_plots = self.study_dir / "plots"

        for d in [
            self.study_dir,
            self.dir_input,
            self.dir_alignment,
            self.dir_tree,
            self.dir_iqtree,
            self.dir_plots
        ]:
            d.mkdir(parents=True, exist_ok=True)

        self.log(f"Created study folder: {self.study_dir}")

    def sanitize_fasta_headers(self, input_fasta, output_fasta):

        def extract_short_label(header: str) -> str:
            if "[" in header and "]" in header:
                return header.split("[", 1)[1].split("]", 1)[0].strip()
            
            parts = header.split()
            if len(parts) >= 3:
                return f"{parts[1]} {parts[2]}"
    
            return header

        mapping = {}
        with open(input_fasta) as fin, open(output_fasta, "w") as fout:
            i = 1
            for line in fin:
                if line.startswith(">"):
                    clean_id = f"seq{i}"
                    original = line[1:].strip()
                    short_label = extract_short_label(original)
                    mapping[clean_id] = short_label
                    fout.write(f">{clean_id}\n")
                    i += 1
                else:
                    fout.write(line)
        return mapping

    def setup_ui(self):
        self.file_frame = ctk.CTkFrame(self.root)
        self.file_frame.pack(pady=20, padx=20, fill="x")
        
        ctk.CTkLabel(self.file_frame, text="FASTA files:", font=ctk.CTkFont(size=16, weight="bold")).pack(pady=10)
        
        self.file_listbox = tk.Listbox(self.file_frame, height=6)
        self.file_listbox.pack(pady=10, padx=10, fill="x")
        
        btn_frame = ctk.CTkFrame(self.file_frame)
        btn_frame.pack(pady=10)
        
        ctk.CTkButton(btn_frame, text="Add FASTA files", command=self.add_files, width=120).pack(side="left", padx=5)
        ctk.CTkButton(btn_frame, text="Clear list", command=self.clear_files, width=120).pack(side="left", padx=5)
        
        self.action_frame = ctk.CTkFrame(self.root)
        self.action_frame.pack(pady=20, padx=20, fill="x")
        
        ctk.CTkLabel(self.action_frame, text="Actions:", font=ctk.CTkFont(size=16, weight="bold")).pack(pady=10)
        
        btn_align_frame = ctk.CTkFrame(self.action_frame)
        btn_align_frame.pack(pady=10, fill="x")
        ctk.CTkButton(btn_align_frame, text="Align sequences", 
                     command=lambda: self.run_threaded(self.align_sequences), 
                     height=50).pack(pady=5, padx=20, fill="x")
        
        btn_tree_frame = ctk.CTkFrame(self.action_frame)
        btn_tree_frame.pack(pady=10, fill="x")
        ctk.CTkButton(btn_tree_frame, text="Phylogenetic analysis", 
                     command=lambda: self.run_threaded(self.phylogenetic_analysis), 
                     height=50).pack(pady=5, padx=20, fill="x")
        
        self.log_frame = ctk.CTkFrame(self.root)
        self.log_frame.pack(pady=20, padx=20, fill="both", expand=True)
        
        ctk.CTkLabel(self.log_frame, text="Logs:", font=ctk.CTkFont(size=14, weight="bold")).pack(pady=10)
        self.log_text = ctk.CTkTextbox(self.log_frame, height=150)
        self.log_text.pack(pady=10, padx=10, fill="both", expand=True)
    
    def log(self, message):
        self.log_text.insert("end", f"{message}\n")
        self.log_text.see("end")
        self.root.update()
    
    def add_files(self):
        files = filedialog.askopenfilenames(
            title="Select FASTA files",
            filetypes=[("FASTA files", "*.fasta *.fas *.fa *.fasta.gz")]
        )
        for file in files:
            if file not in self.fasta_files:
                self.fasta_files.append(file)
                self.file_listbox.insert("end", os.path.basename(file))
        self.log(f"Added {len(files)} files. Total: {len(self.fasta_files)}")
    
    def clear_files(self):
        self.fasta_files.clear()
        self.file_listbox.delete(0, "end")
        self.log("File list cleared")
    
    def run_threaded(self, func):
        thread = threading.Thread(target=func, daemon=True)
        thread.start()
    
    
    def combine_fasta(self, input_files, output_file):
        with open(output_file, "w") as out:
            for fasta_file in input_files:
                with open(fasta_file, "r") as f:
                    out.write(f.read())
        return output_file
    
    def win_path_to_wsl(self, p: str) -> str:
        p = p.replace("\\", "/")
        # C:/Users/... -> /mnt/c/Users/...
        if len(p) >= 3 and p[1:3] == ":/":
            drive = p[0].lower()
            p = "/mnt/" + drive + p[2:]
        return p
    
    def run_mafft(self, input_fasta, output_alignment):

        wsl_input = self.win_path_to_wsl(str(input_fasta))
        # cmd = ["wsl", "mafft", "--quiet", "--thread", "4", "--auto", wsl_input]
        cmd = ["wsl", "mafft", "--quiet", "--thread", "2", "--auto", wsl_input]

        self.log(f"[MAFFT] Running: {' '.join(cmd)}")
        
        with open(output_alignment, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, cwd=self.study_dir)
        
        if proc.returncode != 0:
            self.log(f"[MAFFT] Error: {proc.stderr}")
            raise RuntimeError("MAFFT failed")
        
        self.log(f"[MAFFT] Alignment saved: {output_alignment}")
        return output_alignment
    
    def align_sequences(self):
        
        if not self.fasta_files:
            messagebox.showwarning("Warning", "Add FASTA files first!")
            return
        
        try:
            self.create_study_dirs()
            self.log("Starting alignment...")
            
            raw_fasta = self.dir_input / "combined_raw.fasta"
            self.combine_fasta(self.fasta_files, raw_fasta)

            clean_fasta = self.dir_input / "combined_clean.fasta"
            self.header_map = self.sanitize_fasta_headers(raw_fasta, clean_fasta)

            combined_fasta = clean_fasta

            self.log(f"Combined FASTA: {combined_fasta}")
            
            
            alignment_file = self.dir_alignment / "alignment.fasta"
            self.run_mafft(combined_fasta, alignment_file)


            self.input_for_mafft = combined_fasta
            self.alignment_file = alignment_file



            messagebox.showinfo("Success", f"Alignment saved:\n{alignment_file}")
            
        except Exception as e:
            self.log(f"Error: {str(e)}")
            messagebox.showerror("Error", str(e))
    
    def run_iqtree(self, alignment_file):
        wsl_alignment = self.win_path_to_wsl(str(alignment_file))
        wsl_prefix = self.win_path_to_wsl(str(self.dir_iqtree / "tree_analysis"))

        cmd = [
            "wsl", "iqtree2",
            "-s", wsl_alignment,
            "-m", "MFP",
            "-bb", "1000",
            "-alrt", "1000",
            "-nt", "AUTO",
            "-pre", wsl_prefix,
            "-redo"
        ]

        self.log(f"[IQ-TREE] Running (WSL): {' '.join(cmd)}")

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        if proc.returncode != 0:
            self.log(f"[IQ-TREE] Error:\n{proc.stderr}")
            raise RuntimeError("IQ-TREE failed in WSL")

        self.log("IQ-TREE completed!")
        return Path(str(self.dir_iqtree / "tree_analysis"))

    
    def plot_tree(self, treefile):
        png_file = self.dir_plots / "phylogenetic_tree.png"
        self.log(f"Creating tree plot: {png_file}")
        
        tree = Phylo.read(treefile, "newick")

        plt.figure(figsize=(14, 10))
        Phylo.draw(tree, do_show=False)
        
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        self.log("Tree plot saved!")
        return png_file
    
    def open_file(self, path: Path):
        if platform.system() == "Windows":
            os.startfile(path)
        elif platform.system() == "Darwin":
            subprocess.run(["open", path])
        else:
            subprocess.run(["xdg-open", path])

    def build_label_map(self, fasta_path):
        mapping = {}
        with open(fasta_path) as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip()
                    parts = header.split(" ", 1)
                    acc = parts[0]
                    rest = parts[1] if len(parts) > 1 else ""
                    if "[" in rest and "]" in rest:
                        species = rest.split("[", 1)[1].split("]", 1)[0]
                    else:
                        species = rest if rest else acc
                    mapping[acc] = species
        return mapping

    def relabel_tree(self, treefile, label_map, out_newick):
        tree = Phylo.read(treefile, "newick")
        for clade in tree.get_terminals():
            if clade.name in label_map:
                clade.name = label_map[clade.name]
        Phylo.write(tree, out_newick, "newick")
        return out_newick

    def phylogenetic_analysis(self):
        
        alignment_file = self.alignment_file

        if not alignment_file.exists():
            messagebox.showwarning("Warning", "Run alignment first!")
            return
        
        try:
            self.log("Starting phylogenetic analysis...")
            
            prefix = self.run_iqtree(alignment_file)
            
            treefile = None
            for ext in [".contree", ".treefile"]:
                candidate = Path(str(prefix) + ext)
                if candidate.exists():
                    treefile = candidate
                    break
            
            if not treefile:
                raise FileNotFoundError("Tree file not found")
            
            self.log(f"Tree file: {treefile}")
            
            if not self.input_for_mafft:
                raise RuntimeError("No input FASTA remembered for MAFFT")
            label_map = self.header_map


            labeled_treefile = self.dir_tree / "tree_labeled.treefile"
            self.relabel_tree(treefile, label_map, labeled_treefile)
            self.log(f"Labeled tree file: {labeled_treefile}")

            png = self.plot_tree(labeled_treefile)
            self.open_file(png)

            
            messagebox.showinfo("Success", 
                f"Analysis complete!\n\n"
                f"Output folder: {self.study_dir}\n"
                f"Tree: {treefile}\n"
                f"Plot: {self.dir_plots / 'phylogenetic_tree.png'}")
            
        except Exception as e:
            self.log(f"Error: {str(e)}")
            messagebox.showerror("Error", str(e))
    
    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = PhylogenyApp()
    app.run()
