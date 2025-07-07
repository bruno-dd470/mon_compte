# -*- coding: utf-8 -*-
import os
import sys
import logging
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez
sys.path.insert(0, r"C:\Users\Bruno_Win\Downloads\Spooky_project_20250625")
sys.path.append(os.path.dirname(os.path.dirname(__file__)))  # Add parent directory to path
from spooky_processor import process_tsv, load_ncbi_config, save_ncbi_config

class SpookyProcessorGUI:
    def __init__(self, master=None):
        if master is None:
            import tkinter as tk
            master = tk.Tk()  # Crée une fenêtre principale si aucun master n’est passé
        self.master = master
        self.setup_ui()
        self.setup_logging()
        self.show_readme()

    def setup_ui(self):
        self.master.title("Spooky Processor v2.0")
        self.master.geometry("900x700")
        self.master.resizable(True, True)

        main_frame = ttk.Frame(self.master, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Input Section
        input_frame = ttk.LabelFrame(main_frame, text="Input", padding="10")
        input_frame.pack(fill=tk.X, pady=5)
        ttk.Label(input_frame, text="TSV File:").grid(row=0, column=0, sticky=tk.W)
        self.input_entry = ttk.Entry(input_frame, width=60)
        self.input_entry.grid(row=0, column=1, padx=5)
        ttk.Button(input_frame, text="Browse", command=self.browse_input).grid(row=0, column=2)

        # Output Section
        output_frame = ttk.LabelFrame(main_frame, text="Output", padding="10")
        output_frame.pack(fill=tk.X, pady=5)
        ttk.Label(output_frame, text="Directory:").grid(row=0, column=0, sticky=tk.W)
        self.output_entry = ttk.Entry(output_frame, width=60)
        self.output_entry.grid(row=0, column=1, padx=5)
        ttk.Button(output_frame, text="Browse", command=self.browse_output).grid(row=0, column=2)

        # Options
        options_frame = ttk.LabelFrame(main_frame, text="Options", padding="10")
        options_frame.pack(fill=tk.X, pady=5)
        self.include_gca = tk.BooleanVar(value=False)
        self.include_resonance = tk.BooleanVar(value=True)
        self.max_workers = tk.IntVar(value=4)
        ttk.Checkbutton(options_frame, text="Include GCA_", variable=self.include_gca).grid(row=0, column=0, sticky=tk.W)
        ttk.Checkbutton(options_frame, text="Include Resonance Data", variable=self.include_resonance).grid(row=0, column=1, sticky=tk.W)
        ttk.Label(options_frame, text="Threads:").grid(row=1, column=0, sticky=tk.W)
        ttk.Spinbox(options_frame, from_=1, to=16, textvariable=self.max_workers, width=4).grid(row=1, column=1, sticky=tk.W)

        # Log Area
        log_frame = ttk.LabelFrame(main_frame, text="Processing Log", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.log_text = tk.Text(log_frame, height=15, state=tk.DISABLED)
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Progress & Controls
        control_frame = ttk.Frame(main_frame)
        control_frame.pack(fill=tk.X, pady=5)
        self.progress = ttk.Progressbar(control_frame, mode='indeterminate')
        self.progress.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
        ttk.Button(control_frame, text="Start", command=self.start_processing).pack(side=tk.RIGHT)

    def setup_logging(self):
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        # Handler fichier
        log_path = os.path.join(os.path.dirname(__file__), "spooky_processor.log")
        file_handler = logging.FileHandler(log_path, mode='w', encoding='utf-8')
        file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

        # Handler GUI
        class GuiLogHandler(logging.Handler):
            def __init__(self, text_widget):
                super().__init__()
                self.text_widget = text_widget
            def emit(self, record):
                msg = self.format(record)
                self.text_widget.config(state=tk.NORMAL)
                self.text_widget.insert(tk.END, msg + "\n")
                self.text_widget.see(tk.END)
                self.text_widget.config(state=tk.DISABLED)
        self.gui_handler = GuiLogHandler(self.log_text)
        self.gui_handler.setLevel(logging.INFO)
        self.gui_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(self.gui_handler)

    def browse_input(self):
        path = filedialog.askopenfilename(filetypes=[("TSV Files", "*.tsv")])
        if path:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, path)

    def browse_output(self):
        path = filedialog.askdirectory()
        if path:
            self.output_entry.delete(0, tk.END)
            self.output_entry.insert(0, path)

    def get_credentials(self):
        config = load_ncbi_config()
        if config.get("email"):
            if messagebox.askyesno("Credentials", "Use saved credentials?"):
                return config["email"], config.get("api_key")
        email = simpledialog.askstring("NCBI Email", "Enter your NCBI email:")
        if not email:
            raise ValueError("Email is required for NCBI access")
        api_key = simpledialog.askstring("API Key", "Optional NCBI API key:", show='*')
        save_ncbi_config(email, api_key)
        return email, api_key

    def start_processing(self):
        input_file = self.input_entry.get().strip()
        output_dir = self.output_entry.get().strip()
        if not all([input_file, output_dir]):
            messagebox.showerror("Error", "Please specify both input and output paths.")
            return

        try:
            email, api_key = self.get_credentials()
            Entrez.email = email
            if api_key:
                Entrez.api_key = api_key

            self.progress.start()
            threading.Thread(
                target=self.run_processing,
                args=(input_file, output_dir),
                daemon=True
            ).start()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to start processing: {str(e)}")

    def run_processing(self, input_file, output_dir):
        try:
            process_tsv(
                input_file,
                output_dir,
                include_gca=self.include_gca.get(),
                max_workers=self.max_workers.get(),
                include_resonance=self.include_resonance.get()
            )
            self.master.after(100, lambda: messagebox.showinfo("Success", "Processing completed!"))
        except Exception as e:
            self.master.after(100, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.master.after(100, self.progress.stop)

    def show_readme(self):
        readme_path = os.path.join(os.path.dirname(__file__), "README.txt")
        if os.path.exists(readme_path):
            try:
                with open(readme_path, "r", encoding="utf-8") as f:
                    content = f.read()
                self.log_text.config(state=tk.NORMAL)
                self.log_text.insert(tk.END, f"README:\n{content}\n\n")
                self.log_text.config(state=tk.DISABLED)
            except Exception as e:
                logging.warning(f"Could not read README: {str(e)}")


# Point d'entrée principal
def main():
    root = tk.Tk()
    app = SpookyProcessorGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()