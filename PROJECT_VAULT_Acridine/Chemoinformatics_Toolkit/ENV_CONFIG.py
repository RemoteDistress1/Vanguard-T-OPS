import os
import shutil
from datetime import datetime

# --- CONFIGURATION ---
# This is where your Vanguard project lives
base_dir = "acridine_project"
sub_dirs = ["scripts", "data", "logs", "visuals"]

# --- THE BUILDER ---
# Logic from Sweigart Ch. 9: os.makedirs()
if not os.path.exists(base_dir):
    os.makedirs(base_dir)
    print(f"[+] Vanguard Base Constructed: ./{base_dir}")
else:
    print(f"[!] Base already exists: ./{base_dir}")

# Build the Sub-Sectors
for folder in sub_dirs:
    # os.path.join handles the slashes for Windows vs Linux automatically
    path = os.path.join(base_dir, folder)
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"   [+] Sector Secured: {folder}")

# --- THE LOGBOOK ---
# Logic: Opening a file in 'a' (append) mode creates it if it doesn't exist
log_path = os.path.join(base_dir, "logs", "status.txt")
time_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

with open(log_path, "a") as log_file:
    log_file.write(f"--- OPERATIONS ACTIVE: {time_stamp} ---\n")
    log_file.write("System: Vanguard Therapeutics Digital Dossier initialized.\n")
    log_file.write("Target: Low-Cost Acridine Photoredox Synthesis.\n")

print(f"[+] Mission Log Updated at: {log_path}")
