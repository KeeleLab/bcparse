# compile.py
"""
Name:      compile.py
Author:     CAG
Version:    2.5
Date:       2026/03/20
Refactored: 2026/05/26
"""

# %% Imports
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


# %% Classes

# --- COMPILE MODE ---
class CompileArgSelectionGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Compiler")

        self.args = {}

        # Main frame that expands to fill window
        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        # Define Tkinter variables
        self.ldist_s_var = tk.BooleanVar(value=False)
        self.ldist_g_var = tk.BooleanVar(value=False)
        self.ccheck_var = tk.BooleanVar(value=True)
        self.collapse_to_parent_var = tk.BooleanVar(value=True)
        self.dist_threshold_var = tk.IntVar(value=1)

        self.setup_ui()

    def setup_ui(self):

        # ----------------------------
        # - Required paramters frame -
        # ----------------------------

        # Required selections frame
        self.req_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.req_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        req_label = tk.Label(self.req_frame, text="Required selections")
        req_label.pack(padx=2, pady=2, anchor="w")

        # Input directory
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, pady=2, expand=True)
        tk.Button(frame, text="Browse Analysis Directory", command=self.browse_xlsx_dir).pack(side=tk.LEFT, padx=5)
        self.xlsx_path_label = tk.Label(frame, text="Analysis files directory: Not selected")
        self.xlsx_path_label.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # Output directory
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, pady=2, expand=True)
        tk.Button(frame, text="Browse Output Directory", command=self.browse_out_dir).pack(side=tk.LEFT, padx=5)
        self.out_path_label = tk.Label(frame, text="Output directory: Not selected")
        self.out_path_label.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # Output prefix entry
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, pady=2)
        tk.Label(frame, text="Output file prefix:").pack(side=tk.LEFT, padx=5)
        self.prefix_entry = tk.Entry(frame)
        self.prefix_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)


        # ---------------------------
        # - Optional base csv frame -
        # ---------------------------

        # Optional selections frame for checkbox
        self.optbase_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.optbase_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        base_label = tk.Label(
            self.optbase_frame, text="Optional: add new analyses to prev. compiled data '...all.csv'"
            )
        base_label.pack(padx=2, pady=2, anchor="w")

        # Optional base CSV
        frame = tk.Frame(self.optbase_frame)
        frame.pack(fill=tk.X, pady=2, expand=True)

        tk.Button(
            frame,
            text="Browse Base CSV",
            command=self.browse_base_csv,
        ).pack(side=tk.LEFT, padx=5)

        self.base_csv_label = tk.Label(
            frame,
            text="Base CSV: Not selected",
        )
        self.base_csv_label.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # -------------------------
        # - Optional checks frame -
        # -------------------------

        # Optional selections frame for checkbox
        self.check_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.check_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        opt_label = tk.Label(self.check_frame, text="Optional analyses")
        opt_label.pack(padx=2, pady=2, anchor="w")

        # Spinbox for dist_threshold
        frame = tk.Frame(self.check_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Distance threshold:").pack(side=tk.LEFT, padx=2)

        self.dist_threshold_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=40,
            textvariable=self.dist_threshold_var,
            width=5,
            )
        self.dist_threshold_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)

        # Checkbox for ldist_samp_lvl
        self.ldist_samp_lvl = tk.Checkbutton(
            self.check_frame, 
            text="Rerun sample-level ldist checks & replace previous (for old runs!)", variable=self.ldist_s_var)
        self.ldist_samp_lvl.pack(anchor="w", padx=5, pady=5)

        # Checkbox for ldist_group_lvl
        self.ldist_group_lvl = tk.Checkbutton(
            self.check_frame, 
            text="Run compiled group-level ldist checks", variable=self.ldist_g_var)
        self.ldist_group_lvl.pack(anchor="w", padx=5, pady=5)
        
        # Checkbox for contam_check
        self.contam_check = tk.Checkbutton(
            self.check_frame,
             text="Run contamination check & get matrices", variable=self.ccheck_var)
        self.contam_check.pack(anchor="w", padx=5, pady=5)       

        self.collapse_to_parent = tk.Checkbutton(
            self.check_frame,
            text="Collapse children to parent within each sample",
            variable=self.collapse_to_parent_var,
        )
        self.collapse_to_parent.pack(anchor="w", padx=5, pady=5)

        # Submit button
        submit_btn = tk.Button(self.main_frame, text="Submit", command=self.on_submit)
        submit_btn.pack(pady=10, anchor="w")

    def browse_xlsx_dir(self):
        directory = filedialog.askdirectory(title="Select Directory Containing Analysis Excel Files", initialdir=os.getcwd())
        if directory:
            self.args["xlsx_path"] = directory
            self.xlsx_path_label.config(text=f"Analysis files directory: {directory}")

    def browse_out_dir(self):
        directory = filedialog.askdirectory(title="Select Output Directory", initialdir=os.getcwd())
        if directory:
            self.args["out_path"] = directory
            self.out_path_label.config(text=f"Output directory: {directory}")
    
    def browse_base_csv(self):
        path = filedialog.askopenfilename(
            title="Select previously compiled all.csv (optional)",
            filetypes=[("CSV files", "*.csv")],
            initialdir=os.getcwd(),
        )
        if path:
            self.args["base_csv_path"] = path
            self.base_csv_label.config(text=f"Base CSV: {path}")


    def on_submit(self):
        
        # prefix
        prefix = self.prefix_entry.get().strip()
        if prefix:
            self.args["out_prefix"] = prefix

        # required fields
        required_fields = ["xlsx_path", "out_path", "out_prefix"]
        missing = [f for f in required_fields if f not in self.args or not self.args[f]]
        if missing:
            messagebox.showerror("Missing fields", f"Please select/enter: {', '.join(missing)}")
            return

        # flags
        self.args["ldist_samp_lvl"]  = bool(self.ldist_s_var.get())
        self.args["ldist_group_lvl"] = bool(self.ldist_g_var.get())
        self.args["contam_check"]    = bool(self.ccheck_var.get())
        self.args["collapse_to_parent"] = bool(self.collapse_to_parent_var.get())

        # distance params
        dt = int(self.dist_threshold_var.get())
        self.args["dist_threshold"] = max(0, dt)  # clamp to >= 0

        self.root.quit()


    def run(self):
        self.root.mainloop()
        try:
            self.root.destroy()
        except tk.TclError:
            pass

        return dict(self.args)


#%% Versions
"""
2.5  - Dropped dist backend param
2.41 - Changed defaults and parse gui messages
2.4  - added optional base csv to compile mode
2.3  - changed to accomodate dist_threshold/backend and dualindex params
2.2  - changed to accomodate out_path/out_prefix changes
2.1  - set compile GUI to start from cwd
2026-05-26 - Refactored into bcparse package structure
"""
