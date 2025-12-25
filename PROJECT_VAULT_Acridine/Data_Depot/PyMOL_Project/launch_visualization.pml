# PyMOL Visualization Script for Acridine Project
# Usage: Open this file in PyMOL to load everything automatically.

# 1. Load the data
load Master_Collection.sdf, Candidates

# 2. Apply Organic Chemistry Styling
hide lines
show sticks
set stick_radius, 0.25
util.cbag  # Color carbon atoms green (standard contrast)

# 3. Highlight the Oxygen atoms (Target Interaction)
select oxygens, elem O
color red, oxygens
show spheres, oxygens
set sphere_scale, 0.3

# 4. View Setup
zoom
bg_color white
set ray_shadows, 0
