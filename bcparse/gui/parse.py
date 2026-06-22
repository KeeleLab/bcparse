# parse.py
"""
Name:      parse.py
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
            "sample_path": ("Select sample fastq file", "Fastq files", "*.*", "sample.fq"),
            "runinfo_path": ("Select Runinfo Excel file", "Excel files", "*.*", "runinfo.xlsx"),
            }

        # Initialize BooleanVars with defaults
        self.dualindex_var = tk.BooleanVar(value=True)
        self.filt_mat_ac_var = tk.BooleanVar(value=True)
        self.mask_var = tk.BooleanVar(value=False)
        self.collapse_ambig_var = tk.BooleanVar(value=False)
        self.legacy_format_var = tk.BooleanVar(value=False)
        self.collapse_to_parent_var = tk.BooleanVar(value=True)
        self.append_spike_ref_var = tk.BooleanVar(value=False)

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
        title, file_desc, file_ext, label = self.file_types[file_type]

        path = filedialog.askopenfilename(
            title=title, filetypes=[(file_desc, file_ext)], initialdir=os.getcwd()
            )

        if path:
            self.args[file_type] = path
            self.labels[file_type].config(text=f"{label}: {path}")

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
        for _, (file_type, (_, _, _, label)) in enumerate(self.file_types.items()):
            frame = tk.Frame(self.req_frame)
            frame.pack(fill=tk.X, expand=True, pady=2)

            tk.Button(
                frame,
                text=f"Browse {label}",
                command=lambda ft=file_type: self.browse_file(ft),
                ).pack(
                    side=tk.LEFT, 
                    padx=5
                    )

            self.labels[file_type] = tk.Label(frame, text=f"{label}: Not selected")
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

        # Spike reference checkbox
        frame = tk.Frame(self.req_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.append_spike_ref_checkbox = tk.Checkbutton(
            frame,
            text="Add Spike",
            variable=self.append_spike_ref_var,
            )
        self.append_spike_ref_checkbox.pack(side=tk.LEFT, padx=2)
        
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

        # Checkbox for ambiguous barcode collapse
        frame = tk.Frame(self.fastq_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.collapse_ambig_checkbox = tk.Checkbutton(
            frame,
            text="Collapse single-N child to parent - mask reqd.",
            variable=self.collapse_ambig_var,
            command=self.toggle_mask_from_collapse_ambig # auto-checks mask!
            )
        self.collapse_ambig_checkbox.pack(side=tk.LEFT, padx=2)

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

        # Checkbox for collapse_to_parent
        frame = tk.Frame(self.dist_frame)
        frame.pack(fill=tk.X, expand=True, pady=2)

        self.collapse_to_parent_checkbox = tk.Checkbutton(
            frame,
            text="Collapse ldist children to parent",
            variable=self.collapse_to_parent_var,
            )
        self.collapse_to_parent_checkbox.pack(side=tk.LEFT, padx=2)

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
    
    def toggle_mask_from_collapse_ambig(self):
        """If ambiguous barcode collapse is checked, ensure mask is also checked."""
        if self.collapse_ambig_var.get() and not self.mask_var.get():
            self.mask_var.set(True)

    def on_submit(self):
        """
        Define action for submit button
        """
        # Check if all required files have been selected
        required_fields = ["sample_path", "runinfo_path", "out_path"]

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
        self.args["collapse_ambig"] = self.collapse_ambig_var.get()
        self.args["mask_qual"] = self.mask_qual_var.get()
        self.args["dist_threshold"] = self.dist_threshold_var.get()
        self.args["filt_mat_ac"] = self.filt_mat_ac_var.get()
        self.args["legacy_format"] = self.legacy_format_var.get()
        self.args["collapse_to_parent"] = self.collapse_to_parent_var.get()
        self.args["append_spike_ref"] = self.append_spike_ref_var.get()

        # All required selections are made; close the window
        self.root.quit()


    def run(self) -> dict:
        """
        Run the GUI and return a dict of args
        """
        self.root.mainloop()  # Wait for user interaction

        try:
            self.root.destroy()  # Safely destroy the root window
        except tk.TclError:
            pass  # Ignore if already destroyed

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
