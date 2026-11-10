# screen.py
"""
Name:       screen.py
Author:     CAG
Version:    1.0.0
Date:       2026/11/10

Tkinter input selection for screen mode.
"""

from __future__ import annotations

import os
import tkinter as tk
from tkinter import filedialog, messagebox


class ScreenArgSelectionGUI:
    """Collect the FASTQ, output path, and screen quality settings."""

    def __init__(self) -> None:
        # Create the main window
        self.root = tk.Tk()
        self.root.title("FASTQ Screen")

        # Create a main frame that auto-expands vs root.geometry
        self.main_frame = tk.Frame(self.root)
        self.main_frame.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        # Dictionary to store selections
        self.args = {}

        # Dictionary to store widget labels
        self.labels = {}

        self.out_label = tk.Label(self.root)

        # Define file types with dialog options
        self.file_types = {
            "sample_path": ("Select sample fastq file", "Fastq files", "*.*", "sample.fq"),
        }

        # Initialize BooleanVars with defaults
        self.mask_var = tk.BooleanVar(value=False)

        # Initialize IntVars with defaults
        self.mean_qual_var = tk.IntVar(value=30)
        self.mismatches_var = tk.IntVar(value=1)
        self.dist_threshold_var = tk.IntVar(value=0)
        self.mask_qual_var = tk.IntVar(value=20)

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
        # Open directory dialog and update file path
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

        # -----------------------------
        # - Distance parameters frame -
        # -----------------------------

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
        self.dist_threshold_spinbox.pack(side=tk.LEFT, padx=2)

        # Submit button
        submit_button = tk.Button(
            self.main_frame, text="Submit", command=self.on_submit
            )
        submit_button.pack(pady=10, anchor="w")

    def on_submit(self):
        """
        Define action for submit button
        """
        # Check if all required files have been selected
        required_fields = ["sample_path", "out_path"]

        # Collect missing fields
        missing_fields = [field for field in required_fields if field not in self.args]

        # Return message if missing
        if missing_fields:
            missing_str = ", ".join(missing_fields)
            messagebox.showerror(
                "Missing Selections", f"Please select the following: {missing_str}"
                )
            return

        # Add screen params to args if all required selections are present
        self.args["mean_qual"] = self.mean_qual_var.get()
        self.args["mismatches"] = self.mismatches_var.get()
        self.args["dist_threshold"] = self.dist_threshold_var.get()
        self.args["mask"] = self.mask_var.get()
        self.args["mask_qual"] = self.mask_qual_var.get()

        # All required selections are made; close the window
        self.root.quit()

    def run(self) -> dict:
        """
        Run the GUI and return a dict of args
        """
        self.root.mainloop()  # Wait for user interaction

        try:
            self.root.destroy()  # Safely destroy
        except tk.TclError:
            pass  # Ignore if already destroyed

        return dict(self.args)


#%% Versions
"""
v1.0.0 20261110
- Added dedicated FASTQ/output selection and quality/distance settings for
  standalone Screen mode.
"""
