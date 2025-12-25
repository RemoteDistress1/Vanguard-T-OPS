
# PyMOL Script to Generate Publication Quality Images
# Usage: Open this file in PyMOL

# 1. Load the Structures
load Candidate_Alpha_79.sdf, ligand
load Target_TopoII.pdb, protein

# 2. Setup the View
hide everything
bg_color white
show cartoon, protein
color gray80, protein
set transparency, 0.4, protein

# 3. Highlight the Ligand (Your Drug)
show sticks, ligand
util.cbag ligand  # Color by atom (green carbon)
show spheres, ligand
set sphere_scale, 0.25, ligand

# 4. Zoom to Binding Site
select active_site, byres protein within 5 of ligand
zoom active_site
show sticks, active_site
util.cbac active_site # Color protein atoms cyan

# 5. Ray Trace & Save Image (High Quality)
set ray_shadows, 0
ray 1200, 1200
png Figure_3_Binding_Mode.png

# 6. Re-orient for a second view (Overview)
zoom protein
ray 1200, 1200
png Figure_4_Protein_Overview.png

print "✅ Images Generated: Figure_3_Binding_Mode.png & Figure_4_Protein_Overview.png"
