# gui.py
"""
Name:       gui.py
Author:     CAG
Version:    2.5
Date:       2026/03/20
"""

# %% Imports
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


# %% Classes

# --- PARSE MODE ---
class ArgSelectionGUI:

    def __init__(self, stock_list: list[str]):
        # Create the main window
        self.root = tk.Tk()
        self.root.title("Parser")

        # Create a main frame that auto-expands vs root.geometry
        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        # Store stock list for dropdown
        self.stock_list = stock_list

        # Dictionary to store selections
        self.args = {}

        # Dictionary to store widget labels
        self.labels = {}

        self.out_label = tk.Label(self.root)
        self.stock_var = tk.StringVar(self.root)

        # Define file types with dialog options
        self.file_types = {
            "sample.fq": ("Select sample fastq file", "Fastq files", "*.*"),
            "runinfo.xlsx": ("Select Runinfo Excel file", "Excel files", "*.*"),
            }

        # Initialize BooleanVars with defaults
        self.dualindex_var = tk.BooleanVar(value=True)
        self.filt_mat_ac_var = tk.BooleanVar(value=True)
        self.mask_var = tk.BooleanVar(value=False)
        self.collapse_var = tk.BooleanVar(value=False)
        self.legacy_format_var = tk.BooleanVar(value=False)

        # Initialize IntVars with defaults
        self.mean_qual_var = tk.IntVar(value=30)
        self.mismatches_var = tk.IntVar(value=1)
        self.mask_qual_var = tk.IntVar(value=20)
        self.dist_threshold_var = tk.IntVar(value=1)

        self.setup_ui()

    def browse_file(self, file_type):
        """
        Define action for file selection buttons
        """
        # Open file dialog and update file path
        title, file_desc, file_ext = self.file_types[file_type]

        path = filedialog.askopenfilename(
            title=title, filetypes=[(file_desc, file_ext)], initialdir=os.getcwd()
            )

        if path:
            self.args[file_type] = path
            self.labels[file_type].config(text=f"{file_type}: {path}")

    def browse_output(self):
        """
        Define action for output directory button
        """
        # Open directory dialog for output selection
        path = filedialog.askdirectory(
            title="Select Output Directory", initialdir=os.getcwd()
            )

        if path:
            self.args["out_path"] = path
            self.out_label.config(text=f"output: {path}")


    def setup_ui(self):
        """
        Build the GUI
        """

        # ------------------
        # - Required input -
        # ------------------ 

        # Create a subframe for required selections
        self.req_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.req_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        # Add a label to the subframe
        req_label = tk.Label(self.req_frame, text="Required selections")
        req_label.pack(padx=2, pady=2, anchor="w")

        # Create selection buttons for files
        for _, (file_type, _) in enumerate(self.file_types.items()):
            frame = tk.Frame(self.req_frame)
            frame.pack(fill=tk.X, expand=True, pady=2)

            tk.Button(
                frame,
                text=f"Browse {file_type}",
                command=lambda ft=file_type: self.browse_file(ft),
                ).pack(
                    side=tk.LEFT, 
                    padx=5
                    )

            self.labels[file_type] = tk.Label(frame, text=f"{file_type}: Not selected")
            self.labels[file_type].pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # Output directory
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Button(
            frame, 
            text="Browse output", 
            command=self.browse_output
            ).pack(
                side=tk.LEFT, padx=5
                )

        self.out_label = tk.Label(frame, text="output: Not selected")
        self.out_label.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # Stock selection drop-down
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Select stock:").pack(side=tk.LEFT, padx=5)

        # Variable to hold the stock selection, default 239M
        self.stock_var.set("239M")

        stock_menu = ttk.Combobox(frame, textvariable=self.stock_var, values=self.stock_list)
        stock_menu.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=False)

        # Dual index checkbox
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.dualindex_checkbox = tk.Checkbutton(
            frame, text="Dual-index", variable=self.dualindex_var
            )
        self.dualindex_checkbox.pack(side=tk.LEFT, padx=2)
        
        # ------------------
        # - Fastq QC frame -
        # ------------------

        # Create a subframe for optional selections
        self.fastq_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.fastq_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        # Add a label to the subframe
        fastq_label = tk.Label(self.fastq_frame, text="Optional - Fastq QC")
        fastq_label.pack(padx=2, pady=2, anchor="w")

        # Checkbox for mask
        frame = tk.Frame(self.fastq_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.mask_checkbox = tk.Checkbutton(
            frame, text="Mask low-q bases with N", variable=self.mask_var
            )
        self.mask_checkbox.pack(side=tk.LEFT, padx=2)

        # Checkbox for collapse
        frame = tk.Frame(self.fastq_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.collapse_checkbox = tk.Checkbutton(
            frame,
            text="Collapse single-N child to parent - mask reqd.",
            variable=self.collapse_var,
            command=self.toggle_mask_from_collapse # auto-checks mask!
            )
        self.collapse_checkbox.pack(side=tk.LEFT, padx=2)

        # Spinbox for mask_qual
        frame = tk.Frame(self.fastq_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Mask PHRED below:").pack(side=tk.LEFT, padx=2)

        self.mask_qual_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=40,
            textvariable=self.mask_qual_var,
            width=5,
            )
        self.mask_qual_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)

        # Spinbox for mean_qual
        frame = tk.Frame(self.fastq_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Drop PHRED below:").pack(side=tk.LEFT, padx=2)

        self.mean_qual_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=40,
            textvariable=self.mean_qual_var,
            width=5,
            )
        self.mean_qual_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)

        # Spinbox for mismatches
        frame = tk.Frame(self.fastq_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        tk.Label(frame, text="Allow mismatches?").pack(side=tk.LEFT, padx=2)

        self.mismatches_spinbox = tk.Spinbox(
            frame,
            from_=0,
            to=5,
            textvariable=self.mismatches_var,
            width=5,
            )
        self.mismatches_spinbox.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=False)


        # ----------------------------
        # - Distance paramters frame -
        # ----------------------------

        # Create a subframe for optional selections
        self.dist_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.dist_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        # Add a label to the subframe
        dist_label = tk.Label(self.dist_frame, text="Optional - Distance parameters")
        dist_label.pack(padx=2, pady=2, anchor="w")

        # Spinbox for dist_threshold
        frame = tk.Frame(self.dist_frame)
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

        # --------------------------
        # - Excel formatting frame -
        # --------------------------

        # Create a subframe for optional selections
        self.fmt_frame = tk.Frame(self.main_frame, borderwidth=2, relief="groove")
        self.fmt_frame.pack(side=tk.TOP, fill=tk.X, expand=False, padx=5, pady=5)

        # Add a label to the subframe
        fmt_label = tk.Label(self.fmt_frame, text="Optional - Excel format")
        fmt_label.pack(padx=2, pady=2, anchor="w")

        # Checkbox for legacy_format
        frame = tk.Frame(self.fmt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.legacy_format_checkbox = tk.Checkbutton(
            frame,
            text="Use legacy format",
            variable=self.legacy_format_var,
            )
        self.legacy_format_checkbox.pack(side=tk.LEFT, padx=2)

        # Checkbox for filt_mat_ac
        frame = tk.Frame(self.fmt_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.filt_mat_ac_checkbox = tk.Checkbutton(
            frame,
            text="Filter matrix above-cutoff",
            variable=self.filt_mat_ac_var,
            )
        self.filt_mat_ac_checkbox.pack(side=tk.LEFT, padx=2)

        # Submit button
        submit_button = tk.Button(
            self.main_frame, text="Submit", command=self.on_submit
            )
        submit_button.pack(pady=10, anchor="w")
    
    def toggle_mask_from_collapse(self):
        """If collapse is checked, ensure mask is also checked."""
        if self.collapse_var.get() and not self.mask_var.get():
            self.mask_var.set(True)

    def on_submit(self):
        """
        Define action for submit button
        """
        # Check if all required files have been selected
        required_fields = ["sample.fq", "runinfo.xlsx", "out_path"]

        # Collect missing fields
        missing_fields = [field for field in required_fields if field not in self.args]

        # Return message if missing
        if missing_fields:
            missing_str = ", ".join(missing_fields)
            messagebox.showerror(
                "Missing Selections", f"Please select the following: {missing_str}"
                )
            return

        # Add stock selection & other params to args if all are present
        self.args["stock"] = self.stock_var.get()
        self.args["dualindex"] = self.dualindex_var.get()
        self.args["mean_qual"] = self.mean_qual_var.get()
        self.args["mismatches"] = self.mismatches_var.get()
        self.args["mask"] = self.mask_var.get()
        self.args["collapse"] = self.collapse_var.get()
        self.args["mask_qual"] = self.mask_qual_var.get()
        self.args["dist_threshold"] = self.dist_threshold_var.get()
        self.args["filt_mat_ac"] = self.filt_mat_ac_var.get()
        self.args["legacy_format"] = self.legacy_format_var.get()

        # All required selections are made; close the window
        self.root.quit()


    def run(self) -> tuple[str, str, str, str, int, bool, int, bool, bool, int, int, bool, bool]:
        """
        Run the GUI and return a list of args
        """
        self.root.mainloop()  # Wait for user interaction

        try:
            self.root.destroy()  # Safely destroy the root window
        except tk.TclError:
            pass  # Ignore if already destroyed

        sample_path = str(self.args.get("sample.fq", ""))
        runinfo_path = str(self.args.get("runinfo.xlsx", ""))
        out_path = str(self.args.get("out_path", ""))
        stock = str(self.args.get("stock", ""))
        dualindex = bool(self.args.get("dualindex", False))
        mean_qual = int(self.args.get("mean_qual", 30))
        mismatches = int(self.args.get("mismatches", 1))
        mask = bool(self.args.get("mask", False))
        collapse = bool(self.args.get("collapse", False))
        mask_qual = int(self.args.get("mask_qual", 20))
        dist_threshold = int(self.args.get("dist_threshold", 1))
        filt_mat_ac = bool(self.args.get("filt_mat_ac", False))
        legacy_format = bool(self.args.get("legacy_format", False))

        return (
            sample_path,
            runinfo_path,
            out_path,
            stock,
            dualindex,
            mean_qual,
            mismatches,
            mask,
            collapse,
            mask_qual,
            dist_threshold,
            filt_mat_ac,
            legacy_format,
            )



# --- COMPILE MODE ---
import tkinter as tk
from tkinter import filedialog, messagebox
import os

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

        xlsx_path     = str(self.args.get("xlsx_path", ""))
        out_path      = str(self.args.get("out_path", ""))
        base_csv_path = self.args.get("base_csv_path", None) # Not str() on purpose to preserve None
        out_prefix    = str(self.args.get("out_prefix", ""))

        ldist_samp    = bool(self.args.get("ldist_samp_lvl", False))
        ldist_group   = bool(self.args.get("ldist_group_lvl", False))
        contam_check  = bool(self.args.get("contam_check", False))

        dist_threshold = int(self.args.get("dist_threshold", 1))

        return (
            xlsx_path,
            out_path,
            base_csv_path,
            out_prefix,
            ldist_samp,
            ldist_group,
            contam_check,
            dist_threshold,
        )



#%% Versions
"""
2.42 - Dropped dist backend param
2.41 - Changed defaults and parse gui messages
2.4 - added optional base csv to compile mode
2.3 - changed to accomodate dist_threshold/backend and dualindex params
2.2 - changed to accomodate out_path/out_prefix changes
2.1 - set compile GUI to start from cwd
"""