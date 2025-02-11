import sys
import re
import os

# Originally By Matthieu Chartier adapted to python by Thomas DesCoteaux
# This program generates the mifView pml files

probesInt = []
ca = []
pseudo = []
protgrid = []
probes = []
probesLab = []
grid = []
pbColors = ["aquamarine", "brightorange", "blue", "red", "limegreen", "lightmagenta"]
gridLab = ["200", "150", "100", "050"]
gridColors = ["br2", "palegreen", "hydrogen", "cyan"]
mifViewFolder = "./"
mifName = ""
mifFolder = ""
mifFile = ""
nbpb = 0
grids = 0
gride = 3
ff = ""
ot = ""

# Read command line arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-m":
        mifFile = sys.argv[i + 1]
    elif sys.argv[i] == "-o":
        mifViewFolder = sys.argv[i + 1]
    elif sys.argv[i] == "-t":
        ot = sys.argv[i + 1]
    elif sys.argv[i] == "-h":
        print("##################\nWelcome to pipeIsoMifView\n##################")
        print("-m         <path to mif file>")
        print("-o         <mifView output directory>")
        print("-t         <output tag name>")
        print("-h         <print help menu>")
        sys.exit()

if mifFile == "":
    print("Error: Missing mif File")
    sys.exit()

# if mifViewFolder is empty, you may need to implement the get_dirs function here.
# mifViewFolder = get_dirs("/Users/matthieuchartier/hive/", "mifView") if mifViewFolder == ""

mifName = os.path.splitext(os.path.basename(mifFile))[0]
mifFolder = os.path.dirname(mifFile)


thinness = 0
probes_size = 6
grid_size = 4
try:
    with open(mifFile, "r") as infile:
        grid = [[] for _ in range(grid_size)]
        probes = [[] for _ in range(probes_size)]
        for line in infile:
            if line.startswith('#protein_grid_distance'):
                numbers = re.findall(r"[-+]?(?:\d*\.*\d+)", line)
                thinness = float(numbers[1]) - float(numbers[0])
            if re.match(r'^#ATOM', line):
                match = re.match(r'#ATOM\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]+)\s+([0-9]+)\s+([a-z0-9]{1})\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([0-9]{1})\s+([0-9]{1})$', line, re.I)
                if match and int(match.group(10)) == 1 and match.group(3) == "CA":
                    ca.append(f"{match.group(1)} {match.group(2)} {match.group(3)} {match.group(4)} {match.group(5)} {match.group(6)} {match.group(7)} {match.group(8)}")
                continue
            elif re.match(r'^#PSEUDO', line):
                match = re.match(r'#PSEUDO\s+([a-z]+)\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)', line, re.I)
                if match:
                    pseudo.append(f"{match.group(1)} {match.group(2)} {match.group(3)} {match.group(4)}")
                continue
            elif re.match(r'^#probe\[([0-9]+)\]\s+([0-9a-z\.-]+)', line, re.I):
                probesLab.append(re.search(r'^#probe\[([0-9]+)\]\s+([0-9a-z\.-]+)', line, re.I).group(2))
                nbpb += 1
                continue
            elif re.match(r'^#zip ([0-9]+)', line, re.I):
                grids = gride = int(re.search(r'^#zip ([0-9]+)', line, re.I).group(1)) if int(re.search(r'^#zip ([0-9]+)', line, re.I).group(1)) != -1 else grids
                continue
            elif re.match(r'^#ff ([0-9a-z]+)', line, re.I):
                ff = re.search(r'^#ff ([0-9a-z]+)', line, re.I).group(1)
                continue
            elif re.match(r'^#PG\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)', line):
                protgrid.append(re.search(r'^#PG\s+([\.0-9-]+)\s+([\.0-9-]+)\s+([\.0-9-]+)', line).groups())
                continue
            elif re.match(r'^#', line) or line.strip() == "":
                continue

            line = line.strip()
            info_list = line.split()
            for item_counter, item in enumerate(info_list):
                info_list[item_counter] = float(item)
            # Store vrtx potential interaction
            for i in range(3, len(info_list) - 5, 3):
                pbid = (i // 3) - 1
                g0 = len(info_list) - 5
                g1 = len(info_list) - 4
                g2 = len(info_list) - 3
                g3 = len(info_list) - 2
                probe_items = (info_list[0], info_list[1], info_list[2], info_list[g0], info_list[g1], info_list[g2], info_list[g3], info_list[i + 1], info_list[i + 2])
                if int(info_list[i]) == 1:
                    for item in probe_items:
                        probes[pbid].append(item)

            # Store vrtx grid presence
            for i in range(len(info_list) - 5, len(info_list)-1):
                gid = i - (len(info_list) - 5)
                gp = info_list[i]
                bu = info_list[-1]
                if int(gp) == 1:
                    grid[gid].append(info_list[0])
                    grid[gid].append(info_list[1])
                    grid[gid].append(info_list[2])
                    grid[gid].append(bu)
except IOError:
    print(f"Can't open mif file {mifFile}")
    sys.exit()

with open(os.path.join(mifViewFolder, f"{mifName}{ot}.pml"), "w") as npml:
    # Print protein
    npml.write("feedback disable,all,output\n")
    with open(os.path.join(mifFolder, f"{mifName}_cpy.pdb"), "r") as pdb_file:
        npml.write('cmd.read_pdbstr("""')
        for line in pdb_file:
            line = line.strip()
            npml.write(line + "\\\n")
        npml.write("TER \\\n")
        npml.write(f'""","{mifName}")\n')
        npml.write('set connect_mode,1\n')

    # Print pseudocenters
    it = 0
    if pseudo:
        npml.write('cmd.read_pdbstr("""')
        for p in pseudo:
            s = p.split()
            npml.write(f"HETATM{it:5d}  N   {s[0].upper():3s} A0000    {float(s[1]):8.3f}{float(s[2]):8.3f}{float(s[3]):8.3f}  0.00 10.00           N\\\n")
            it += 1
        npml.write("TER \\\n")
        npml.write('""","pseudocenters")\n')
        npml.write("hide nonbonded\n")
        npml.write("set connect_mode,1\n")

    # Print mif points for each probe and each grid resolution
    import numpy as np
    probesInt = np.zeros((nbpb, gride+1, 2))
    for i in range(nbpb):  # Loop each probe
        if len(probes[i]) > 0:  # Check if the probe has elements
            for g in range(grids, gride + 1):  # Loop each grid resolution
                probesInt[i][g][0] = it
                npml.write('cmd.read_pdbstr("""')
                for j in range(0, len(probes[i]), 9):  # For each node
                    if probes[i][j + 3 + g] == 1:  # If it's in this grid resolution
                        npml.write(f'HETATM{it:5d}  N   {i:3} A0000    {probes[i][j]:8.3f}{probes[i][j + 1]:8.3f}{probes[i][j + 2]:8.3f}  0.00 10.00           N\\\n')
                        if it != 99999:
                            it += 1
                npml.write(f'TER \\\n')
                npml.write(f'""","{probesLab[i]}_{gridLab[g]}")\n')
                probesInt[i][g][1] = it - 1

    #print grid points
    for i in range(grids, gride + 1):
        it = 0
        if i < len(grid) and len(grid[i]) > 0:
            npml.write('cmd.read_pdbstr("""')
            for j in range(0, len(grid[i]), 4):
                # Ensure we have enough elements to unpack
                if j + 3 < len(grid[i]):
                    npml.write(f"HETATM{it:5d}  N   {gridLab[i]} A0000    "
                               f"{grid[i][j]:8.3f}{grid[i][j + 1]:8.3f}{grid[i][j + 2]:8.3f}  "
                               f"0.00{grid[i][j + 3]:6.2f}           N\\\n")
                    it += 1
                    if it == 99999:
                        break
            npml.write("TER \\\n")
            npml.write(f'""","{gridLab[i]}")\n')
    npml.write("\n")

    # print protein grid points
    if len(protgrid) > 0:
        npml.write('cmd.read_pdbstr("""\n')
        for i in range(len(protgrid)):
            s = protgrid[i].split()
            npml.write(f"HETATM{0:5d}  N   {'PGD':3s} A0000    "
                       f"{float(s[0]):8.3f}{float(s[1]):8.3f}{float(s[2]):8.3f}  "
                       f"0.00 10.00           N\\\n")
        npml.write("TER \\\n\"\"\", \"" + mifName + "_grid\")\n")

    npml.write("\n")
    npml.write("feedback enable,all,output\norient\nshow cartoon, " + mifName + "\n")
    npml.write("remove (resn HOH)\nshow sticks, HET & " + mifName + "\n")
    npml.write("color white," + mifName + "_grid\n")
    npml.write("show nonbonded," + mifName + "_grid\n")
    npml.write("\n")

    for i in range(nbpb):
        if len(probes[i]) > 0:
            for g in range(grids, gride + 1):
                if probesInt[i][g][0] != probesInt[i][g][1]:
                    npml.write(f"show spheres, {probesLab[i]}_{gridLab[g]}\n")
                    npml.write(f"set sphere_scale,0.2,{probesLab[i]}_{gridLab[g]}\n")
                    npml.write("rebuild\n")
                    npml.write(f"color {pbColors[i]},{probesLab[i]}_{gridLab[g]}\n")
                    npml.write(f"hide nonbonded,{probesLab[i]}_{gridLab[g]}\n\n")

    npml.write("\n")

    for i in range(3):
        if i < len(grid) and len(grid[i]) > 0:
            if grid[i][0] != grid[i][1]:
                npml.write(f"color {gridColors[i]},{gridLab[i]}\n")
                npml.write(f"show nonbonded,{gridLab[i]}\n")

    # Finalize
    npml.write("\n")
    npml.write("set sphere_scale, 0.3, pseudocenter\n")
    npml.write("color aquamarine, resn HYD & pseudocenters\n")
    npml.write("color brightorange, resn ARM & pseudocenters\n")
    npml.write("color blue, resn DON & pseudocenters\n")
    npml.write("color red, resn ACC & pseudocenters\n")
    npml.write("color limegreen, resn DOA & pseudocenters\n")
    npml.write("show spheres, pseudocenters\n")

print("File has been generated!")